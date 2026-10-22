"""Testable ABC Atlas adapter; scientific logic is independent of R and network.

The upstream SDK owns manifests, file downloads and hash validation. This module
owns explicit table joins, exact identifiers, H5AD subsetting and provenance.
No expression file is downloaded without allow_downloads=True. No network calls
or environment installation happen at import time. See docs/ARCHITECTURE.md.
"""
from __future__ import annotations

import contextlib
import errno
import hashlib
import importlib.metadata
import json
import math
import numbers
import os
from pathlib import Path
import platform
import re
import shutil
import sys
import time
import traceback
import warnings

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

VERSION = '0.2.0.9007'
DEFAULT_RELEASE = '20260711'
REGISTRY = json.loads(Path(__file__).with_name('registry.json').read_text())
KINDS = {'metadata': 'list_metadata_files',
         'expression_matrices': 'list_expression_matrix_files',
         'image_volumes': 'list_image_volume_files',
         'mapmycells': 'list_mapmycells_files'}
KEY_PATTERN = re.compile(r'(^cell_label$|label$|identifier$|alias$|barcode$|_id$|^label$|^exp_component_name$)')


def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def positive_int(x, name):
    if isinstance(x, bool) or not isinstance(x, (int, np.integer)) or x < 1:
        raise ValueError(f'{name} must be a positive integer.')
    return int(x)


def exact_directory(requested, directories):
    """Resolve case only against an existing manifest; never invent a path."""
    if requested in directories:
        return requested
    matches = [d for d in directories if d.casefold() == requested.casefold()]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(f'Directory {requested!r} absent or ambiguous. Use list_datasets(live=TRUE).')


def _nonnegative_number(value, name, minimum=0.0):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be numeric.') from exc
    if not math.isfinite(number) or number < minimum:
        raise ValueError(f'{name} must be finite and >= {minimum}.')
    return number


def _manifest_name(release):
    if release == 'latest':
        return None
    if re.fullmatch(r'\d{8}', str(release)):
        return f'releases/{release}/manifest.json'
    if re.fullmatch(r'releases/\d{8}/manifest\.json', str(release)):
        return str(release)
    raise ValueError('release must be YYYYMMDD, a releases/YYYYMMDD/manifest.json path, or latest.')


def _network_policy(retries=5, backoff_seconds=2.0, connect_timeout_seconds=30,
                    read_timeout_seconds=300, chunk_mb=8):
    retries = int(retries)
    if retries < 0:
        raise ValueError('network_retries must be >= 0.')
    return {
        'network_retries': retries,
        'network_backoff_seconds': _nonnegative_number(backoff_seconds, 'network_backoff_seconds'),
        'network_connect_timeout_seconds': int(_nonnegative_number(
            connect_timeout_seconds, 'network_connect_timeout_seconds', minimum=1)),
        'network_read_timeout_seconds': int(_nonnegative_number(
            read_timeout_seconds, 'network_read_timeout_seconds', minimum=1)),
        'download_chunk_bytes': int(math.ceil(_nonnegative_number(
            chunk_mb, 'download_chunk_mb', minimum=0.001) * 1024**2)),
    }


def _is_transient_network_error(error):
    if isinstance(error, (BrokenPipeError, ConnectionError, TimeoutError)):
        return True
    name = type(error).__name__.casefold()
    message = str(error).casefold()
    transient_names = ('endpointconnection', 'connectionclosed', 'readtimeout',
                       'connecttimeout', 'incompleteread', 'connectionreset',
                       'responsestreaming')
    transient_text = ('could not connect to the endpoint', 'connection was closed',
                      'read timeout', 'timed out', 'broken pipe', 'response stream',
                      'connection reset', 'temporarily unavailable', 'slowdown',
                      'service unavailable', 'incomplete read', 'request timeout',
                      'internal error', 'bad gateway', 'gateway timeout',
                      'too many requests')
    if any(token in name for token in transient_names) or any(
            token in message for token in transient_text):
        return True
    # Botocore ClientError exposes structured status/code information. Retry
    # throttling, server and timeout responses, but never hide permanent 4xx
    # errors such as AccessDenied or NoSuchKey.
    response = getattr(error, 'response', None)
    if isinstance(response, dict):
        metadata = response.get('ResponseMetadata', {}) or {}
        status = metadata.get('HTTPStatusCode')
        code = str((response.get('Error', {}) or {}).get('Code', '')).casefold()
        if status in (408, 429, 500, 502, 503, 504):
            return True
        if code in {'requesttimeout', 'requesttimeoutexception', 'slowdown',
                    'throttling', 'throttlingexception', 'internalerror',
                    'serviceunavailable'}:
            return True
    return False


def _retry_call(operation, policy, context):
    attempts = int(policy.get('network_retries', 0)) + 1
    delay = float(policy.get('network_backoff_seconds', 0))
    last = None
    for attempt in range(1, attempts + 1):
        try:
            return operation()
        except Exception as error:
            last = error
            if not _is_transient_network_error(error) or attempt >= attempts:
                raise
            wait = min(delay * (2 ** (attempt - 1)), 60.0)
            print(f'importABCatlas: transient network failure during {context}; '
                  f'retrying {attempt}/{attempts - 1} in {wait:.1f}s: {error}',
                  file=sys.stderr, flush=True)
            if wait:
                time.sleep(wait)
    raise last


def _resilient_s3_cache_class(pinned_manifest, policy):
    """Create an SDK-compatible cache with pinned-manifest and resumable downloads."""
    try:
        import boto3
        import tqdm
        from botocore import UNSIGNED
        from botocore.config import Config
        from abc_atlas_access.abc_atlas_cache.cloud_cache import S3CloudCache
        from abc_atlas_access.abc_atlas_cache.utils import file_hash_from_path
    except ImportError as error:
        raise ImportError('Install the Python dependencies with setup_importABCatlas().') from error

    class ImportABCS3CloudCache(S3CloudCache):
        def __init__(self, cache_dir, bucket_name, ui_class_name=None, auth_required=False):
            self._importabcatlas_pinned_manifest = pinned_manifest
            self._importabcatlas_network_policy = dict(policy)
            self._s3_client = None
            super().__init__(cache_dir=cache_dir, bucket_name=bucket_name,
                             ui_class_name=ui_class_name, auth_required=auth_required)

        @property
        def s3_client(self):
            if self._s3_client is None:
                kwargs = {
                    'connect_timeout': self._importabcatlas_network_policy[
                        'network_connect_timeout_seconds'],
                    'read_timeout': self._importabcatlas_network_policy[
                        'network_read_timeout_seconds'],
                    'retries': {
                        'max_attempts': self._importabcatlas_network_policy[
                            'network_retries'] + 1,
                        'mode': 'standard',
                    },
                }
                if not self._auth_required:
                    kwargs['signature_version'] = UNSIGNED
                self._s3_client = boto3.client('s3', config=Config(**kwargs))
            return self._s3_client

        def _list_all_manifests(self):
            # An explicit release is an exact S3 object. Avoid a fragile and
            # unnecessary bucket-wide releases/ listing on every cold cache.
            if self._importabcatlas_pinned_manifest is not None:
                return [self._importabcatlas_pinned_manifest]
            return _retry_call(
                lambda: S3CloudCache._list_all_manifests(self),
                self._importabcatlas_network_policy,
                'listing available manifests')

        def _object_size(self, key):
            response = _retry_call(
                lambda: self.s3_client.head_object(Bucket=self._bucket_name, Key=str(key)),
                self._importabcatlas_network_policy,
                f'checking exact object size for {key}')
            return int(response['ContentLength'])

        def _stream_object(self, key, destination, expected_hash=None):
            """Download one S3 object with exact-size validation and resume support.

            Interrupted bytes are retained as ``.part`` and resumed with an S3
            Range request. A completed object is moved into place only after its
            exact byte count and, when available, manifest hash are verified.
            """
            destination = Path(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            part = Path(str(destination) + '.part')
            expected_size = self._object_size(key)

            def close_response(response):
                if response is None:
                    return
                body = response.get('Body') if isinstance(response, dict) else None
                if body is not None and hasattr(body, 'close'):
                    with contextlib.suppress(Exception):
                        body.close()

            def verify_complete(path):
                if not path.is_file() or path.stat().st_size != expected_size:
                    return False
                if expected_hash is None:
                    return True
                return file_hash_from_path(path) == expected_hash

            # A verified destination needs no transfer. An unverified short file
            # is adopted as a resumable partial; an unverified full/corrupt file
            # is discarded because appending cannot repair it.
            if destination.exists():
                if verify_complete(destination):
                    return destination
                if destination.stat().st_size < expected_size and not part.exists():
                    destination.replace(part)
                else:
                    destination.unlink(missing_ok=True)
            if part.exists() and part.stat().st_size > expected_size:
                part.unlink()

            attempts = self._importabcatlas_network_policy['network_retries'] + 1
            delay = self._importabcatlas_network_policy['network_backoff_seconds']
            chunk_size = self._importabcatlas_network_policy['download_chunk_bytes']
            for attempt in range(1, attempts + 1):
                # A previous attempt may have written the last byte but failed
                # while closing the response. Verify/finalize without requesting
                # an invalid zero-length Range.
                if part.exists() and part.stat().st_size == expected_size:
                    if verify_complete(part):
                        os.replace(part, destination)
                        return destination
                    part.unlink(missing_ok=True)

                offset = part.stat().st_size if part.exists() else 0
                request = {'Bucket': self._bucket_name, 'Key': str(key)}
                if offset:
                    request['Range'] = f'bytes={offset}-'
                response = None
                body = None
                progress = None
                try:
                    response = self.s3_client.get_object(**request)
                    status = int(response.get('ResponseMetadata', {}).get('HTTPStatusCode', 200))
                    content_range = str(response.get('ContentRange', ''))
                    valid_range = (status == 206 and
                                   content_range.casefold().startswith(
                                       f'bytes {offset}-'.casefold()))
                    if offset and not valid_range:
                        # Some S3-compatible endpoints may ignore or return the
                        # wrong Range. Close that response before restarting from
                        # zero; appending it would silently corrupt the object.
                        close_response(response)
                        response = None
                        part.unlink(missing_ok=True)
                        offset = 0
                        response = self.s3_client.get_object(
                            Bucket=self._bucket_name, Key=str(key))
                    body = response.get('Body')
                    if body is None:
                        raise IOError(f'S3 returned no response body for {key}.')
                    progress = tqdm.tqdm(
                        desc=Path(str(key)).name,
                        total=expected_size,
                        initial=offset,
                        unit='B', unit_scale=True, unit_divisor=1024,
                        leave=True)
                    mode = 'ab' if offset else 'wb'
                    with open(part, mode) as stream:
                        for chunk in body.iter_chunks(chunk_size=chunk_size):
                            if chunk:
                                stream.write(chunk)
                                progress.update(len(chunk))
                        stream.flush()
                        with contextlib.suppress(OSError):
                            os.fsync(stream.fileno())
                    actual = part.stat().st_size
                    if actual != expected_size:
                        raise IOError(f'Incomplete response for {key}: received {actual:,} of '
                                      f'{expected_size:,} bytes.')
                    if expected_hash is not None and file_hash_from_path(part) != expected_hash:
                        # A hash mismatch cannot be repaired by resuming. Remove
                        # the corrupted bytes and retry from zero.
                        part.unlink(missing_ok=True)
                        raise IOError(f'Checksum mismatch for {key}.')
                    os.replace(part, destination)
                    return destination
                except Exception as error:
                    retryable = (_is_transient_network_error(error) or
                                 'Incomplete response' in str(error) or
                                 'Checksum mismatch' in str(error))
                    if not retryable:
                        raise
                    if attempt >= attempts:
                        raise RuntimeError(
                            f'Could not download {key} after {attempts} attempts. '
                            f'A resumable partial file is retained at {part}. Last error: {error}') from error
                    wait = min(delay * (2 ** (attempt - 1)), 60.0)
                    print(f'importABCatlas: interrupted download of {key}; retained '
                          f'{part.stat().st_size if part.exists() else 0:,} bytes and will '
                          f'resume in {wait:.1f}s ({attempt}/{attempts - 1}).',
                          file=sys.stderr, flush=True)
                    if wait:
                        time.sleep(wait)
                finally:
                    if progress is not None:
                        progress.close()
                    close_response(response)
            raise RuntimeError(f'Could not download {key}.')

        def _download_manifest(self, manifest_name):
            path = self._cache_dir / manifest_name
            self._stream_object(manifest_name, path, expected_hash=None)

        def _download_file(self, file_attributes, force_download=False,
                           skip_hash_check=False):
            verified = (file_attributes.local_path.exists() and
                        self._check_successful_download(file_attributes))
            if verified and not force_download:
                return False
            if force_download:
                # The upstream wrapper also sets force_download for an existing
                # but unverified partial file. Remove only its success record;
                # _stream_object decides whether the bytes are resumable.
                self._remove_download_from_list(file_attributes)
            expected_hash = None if skip_hash_check else file_attributes.file_hash
            self._stream_object(file_attributes.relative_path,
                                file_attributes.local_path,
                                expected_hash=expected_hash)
            self._update_list_of_downloads(file_attributes)
            return True

    return ImportABCS3CloudCache


def create_cache(download_base, release=DEFAULT_RELEASE, offline=False,
                 check_disk_space=True, disk_reserve_gb=1.0,
                 disk_safety_factor=1.05, allow_unknown_size=False,
                 network_retries=5, network_backoff_seconds=2.0,
                 network_connect_timeout_seconds=30,
                 network_read_timeout_seconds=300, download_chunk_mb=8):
    try:
        from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
    except ImportError as e:
        raise ImportError('Install the Python dependencies with setup_importABCatlas().') from e
    base = Path(download_base).expanduser().resolve()
    manifest_name = _manifest_name(release)
    if offline and not base.exists():
        raise FileNotFoundError(f'Offline cache does not exist: {base}')
    if offline and manifest_name is None:
        raise ValueError('Offline mode requires an explicit release, not latest.')
    base.mkdir(parents=True, exist_ok=True)
    policy = _network_policy(network_retries, network_backoff_seconds,
                             network_connect_timeout_seconds,
                             network_read_timeout_seconds, download_chunk_mb)
    if offline:
        cache = AbcProjectCache.from_local_cache(base)
        cache.load_manifest(manifest_name)
    else:
        cache_class = _resilient_s3_cache_class(manifest_name, policy)
        lower = cache_class(base, bucket_name=AbcProjectCache._default_bucket_name,
                            auth_required=False, ui_class_name='AbcProjectCache')
        cache = _retry_call(lambda: AbcProjectCache(lower), policy,
                            'loading the Allen manifest')
        if manifest_name is None:
            cache.load_latest_manifest()
    cache._local = bool(offline)
    factor = _nonnegative_number(disk_safety_factor, 'disk_safety_factor', minimum=1.0)
    reserve = _nonnegative_number(disk_reserve_gb, 'disk_reserve_gb')
    cache._importabcatlas_download_base = base
    cache._importabcatlas_network_policy = policy
    cache._importabcatlas_disk_policy = {
        'check_disk_space': bool(check_disk_space),
        'disk_reserve_bytes': int(math.ceil(reserve * 1024**3)),
        'disk_reserve_gb': reserve,
        'disk_safety_factor': factor,
        'allow_unknown_size': bool(allow_unknown_size),
    }
    return cache


def files(cache, directory, kind='metadata'):
    fn = getattr(cache, KINDS[kind], None)
    if fn is None:
        raise RuntimeError(f'Upstream SDK lacks {KINDS[kind]}; update abc_atlas_access.')
    # Allen's SDK raises when an existing directory has an empty optional
    # sub-directory. Treat only that documented/observed condition as empty;
    # propagate unrelated ValueError/KeyError failures so network/manifest bugs
    # are never hidden.
    try:
        return list(fn(directory))
    except Exception as error:
        # abc_atlas_access has raised different exception classes for this same
        # legitimate condition across releases. Suppress by exact message shape,
        # never by exception class, so network and malformed-manifest errors pass through.
        message = str(error).lower()
        empty_category = (
            'sub-directory is empty' in message or
            ('no ' in message and ' files found in directory ' in message)
        )
        if empty_category:
            return []
        raise


def resolve_spec(cache, dataset, overrides=None):
    overrides = overrides or {}
    if dataset in REGISTRY:
        spec = dict(REGISTRY[dataset])
    else:
        spec = dict(directory=dataset,cell_file='cell_metadata',gene_directory=dataset,
                    gene_file='gene',taxonomy=None,spatial=False,expression_directories=[dataset])
    allowed = {'directory','cell_file','gene_directory','gene_file','taxonomy',
               'expression_directories','spatial','file_selector','imputed','mapping_directory',
               'default_data_type','available_data_types','modality'}
    if set(overrides)-allowed:
        raise ValueError(f'Unknown override fields: {sorted(set(overrides)-allowed)}')
    spec.update(overrides)
    directories = list(cache.list_directories)
    for name in ('directory','gene_directory','taxonomy','mapping_directory'):
        if spec.get(name):
            spec[name] = exact_directory(spec[name], directories)
    spec['expression_directories'] = [exact_directory(x,directories)
        for x in as_list(spec['expression_directories']) if any(d.casefold()==x.casefold() for d in directories)]
    if not spec['expression_directories']:
        spec['expression_directories'] = [spec['directory']]
    selector = spec.get('file_selector')
    if selector:
        # Native HMBA species share one cell_metadata table.  Species-specific
        # filtering is performed with donor.species_genus before previewing
        # cells.  Gene tables, however, are exact species-specific assets.
        available_cells = files(cache, spec['directory'])
        if 'cell_file' not in overrides:
            if spec['cell_file'] not in available_cells:
                candidates = [x for x in available_cells
                              if selector.casefold() in x.casefold()
                              and 'cell_metadata' in x.casefold()]
                if len(candidates) != 1:
                    raise ValueError(
                        f'{dataset}: cannot identify shared cell metadata; '
                        f'available: {available_cells}. Supply cell_file override only after '
                        'checking the manifest.'
                    )
                spec['cell_file'] = candidates[0]
        if 'gene_file' not in overrides:
            available_genes = files(cache, spec['gene_directory'])
            expected_gene_file = f'{selector.casefold()}_gene'
            exact = [x for x in available_genes
                     if x.casefold() == expected_gene_file.casefold()]
            if len(exact) != 1:
                raise ValueError(
                    f'{dataset}: expected exactly one species gene table '
                    f'{expected_gene_file!r}; available: {available_genes}. '
                    'Example-expression tables are not valid gene metadata. '
                    'Use gene_file override only after checking the manifest.'
                )
            spec['gene_file'] = exact[0]
        spec['species_filter'] = selector.capitalize()
    return spec


def _first_attribute(obj, names):
    if obj is None:
        return None
    for name in names:
        if isinstance(obj, dict) and name in obj:
            value = obj[name]
        else:
            value = getattr(obj, name, None)
        if value is not None:
            return value
    return None


def _non_downloading_file_record(cache, directory, file_name):
    """Inspect a manifest/local path without asking AbcProjectCache to download."""
    inner = getattr(cache, 'cache', None)
    if inner is not None and hasattr(inner, 'get_file_path'):
        try:
            record = inner.get_file_path(directory=directory, file_name=file_name)
        except TypeError:
            record = inner.get_file_path(directory, file_name)
        return record, inner
    # Synthetic fixtures and direct local-cache adapters are safe to query.
    if getattr(cache, '_local', False) and hasattr(cache, 'get_file_path'):
        try:
            record = cache.get_file_path(directory=directory, file_name=file_name)
        except TypeError:
            record = cache.get_file_path(directory, file_name)
        return record, cache
    return None, inner


def _asset_file_info(cache, directory, file_name):
    """Return exact per-file size/cache information without downloading.

    The underlying Allen cache exposes a non-downloading ``get_file_path``
    record containing the intended local path and manifest file attributes.
    For uncached S3 files, exact ``ContentLength`` is obtained for the selected
    object only.  Directory totals are never substituted for file size.
    """
    result = {
        'file_bytes': None,
        'size_source': 'unknown',
        'cached': False,
        'download_required_bytes': None,
        'partial_bytes': 0,
        'relative_path': None,
        'size_error': None,
    }
    try:
        record, inner = _non_downloading_file_record(cache, directory, file_name)
    except Exception as error:
        result['size_error'] = f'cache inspection failed: {error}'
        return result

    attrs = None
    local_path = None
    exists = False
    if isinstance(record, dict):
        attrs = record.get('file_attributes')
        local_path = record.get('local_path')
        exists = bool(record.get('exists', False))
    elif record is not None:
        local_path = record

    if attrs is not None:
        relative = _first_attribute(attrs, ('relative_path', 'remote_path', 'key'))
        if relative is not None:
            result['relative_path'] = str(relative)
        size = _first_attribute(attrs, ('file_size', 'size_bytes', 'content_length',
                                        'ContentLength', 'file_size_bytes', 'size'))
        try:
            numeric_size = int(size) if size is not None and not isinstance(size, bool) else None
        except (TypeError, ValueError, OverflowError):
            numeric_size = None
        if numeric_size is not None and numeric_size >= 0:
            result['file_bytes'] = numeric_size
            result['size_source'] = 'manifest_attributes'

    local_size = None
    if local_path is not None:
        local_path = Path(local_path).expanduser()
        exists = exists or local_path.is_file()
        if local_path.is_file():
            local_size = int(local_path.stat().st_size)

    is_local_cache = bool(getattr(cache, '_local', False))
    cached = bool(exists)
    # S3 cache records may point at an incomplete/unverified local file. Mirror
    # the SDK's integrity decision when that private helper is available.
    if cached and not is_local_cache and attrs is not None:
        checker = getattr(inner, '_check_successful_download', None)
        if callable(checker):
            try:
                cached = bool(checker(attrs))
            except Exception as error:
                cached = False
                result['size_error'] = f'cache verification failed: {error}'
    result['cached'] = cached

    # A verified local file is the most direct exact size. Never use a partial
    # or unverified local file as the expected remote download size.
    if cached and local_size is not None:
        result['file_bytes'] = local_size
        result['size_source'] = 'local_file'

    if not cached and result['relative_path'] and not is_local_cache:
        client = getattr(inner, 's3_client', None) or getattr(inner, '_s3_client', None)
        bucket = (getattr(inner, 'bucket_name', None) or getattr(inner, '_bucket_name', None)
                  or getattr(inner, 'bucket', None))
        key = result['relative_path']
        # Be defensive if an SDK version exposes an s3:// URI rather than a key.
        if key.startswith('s3://'):
            without_scheme = key[5:]
            uri_bucket, _, uri_key = without_scheme.partition('/')
            bucket = uri_bucket or bucket
            key = uri_key
        if client is not None and bucket and key:
            policy = getattr(cache, '_importabcatlas_network_policy', {
                'network_retries': 0, 'network_backoff_seconds': 0,
            })
            try:
                response = _retry_call(
                    lambda: client.head_object(Bucket=bucket, Key=key),
                    policy,
                    f'checking exact file size for {directory}/{file_name}')
                result['file_bytes'] = int(response['ContentLength'])
                result['size_source'] = 's3_head_object'
                result['size_error'] = None
            except Exception as head_error:
                # Some endpoints reject HEAD but permit exact-key LIST. Retry
                # transient failures, but never substitute a directory total.
                try:
                    response = _retry_call(
                        lambda: client.list_objects_v2(Bucket=bucket, Prefix=key),
                        policy,
                        f'listing exact file key for {directory}/{file_name}')
                    exact = [obj for obj in response.get('Contents', [])
                             if obj.get('Key') == key]
                    if len(exact) == 1:
                        result['file_bytes'] = int(exact[0]['Size'])
                        result['size_source'] = 's3_list_objects_v2'
                        result['size_error'] = None
                    elif result['file_bytes'] is None:
                        result['size_error'] = f'cannot obtain exact remote size: {head_error}'
                except Exception as list_error:
                    if result['file_bytes'] is None:
                        result['size_error'] = (f'cannot obtain exact remote size: HEAD={head_error}; '
                                                f'LIST={list_error}')

    partial_bytes = 0
    if not cached and local_path is not None and result['file_bytes'] is not None:
        partial_path = Path(str(local_path) + '.part')
        if partial_path.is_file():
            partial_bytes = min(int(partial_path.stat().st_size), int(result['file_bytes']))
        elif local_size is not None and local_size < int(result['file_bytes']):
            # Adoptable interrupted downloads made by older SDK versions.
            partial_bytes = int(local_size)
    result['partial_bytes'] = partial_bytes
    if cached:
        result['download_required_bytes'] = 0
    elif result['file_bytes'] is not None:
        result['download_required_bytes'] = max(0, int(result['file_bytes']) - partial_bytes)
    return result

def _disk_status(cache, download_bytes):
    policy = getattr(cache, '_importabcatlas_disk_policy', {})
    root = Path(getattr(cache, '_importabcatlas_download_base', Path.cwd()))
    root.mkdir(parents=True, exist_ok=True)
    free = int(shutil.disk_usage(root).free)
    reserve = int(policy.get('disk_reserve_bytes', 0))
    factor = float(policy.get('disk_safety_factor', 1.0))
    required = None if download_bytes is None else int(math.ceil(int(download_bytes) * factor) + reserve)
    return {
        'cache_root': str(root),
        'free_bytes': free,
        'disk_reserve_bytes': reserve,
        'disk_safety_factor': factor,
        'required_free_bytes': required,
        'enough_space': None if required is None else bool(free >= required),
    }


def _assert_download_space(cache, download_bytes, context, unknown_files=None):
    policy = getattr(cache, '_importabcatlas_disk_policy', {})
    if not policy.get('check_disk_space', True) or download_bytes == 0:
        return _disk_status(cache, download_bytes)
    status = _disk_status(cache, download_bytes)
    if download_bytes is None:
        if not policy.get('allow_unknown_size', False):
            files_text = ', '.join(unknown_files or []) or 'one or more files'
            raise RuntimeError(
                f'Cannot verify free disk space before {context}: exact size is unknown for {files_text}. '
                'Retry when S3 file-size lookup is available, or deliberately set '
                'options(importABCatlas.allow_unknown_size=TRUE) to accept the risk.'
            )
        return status
    if not status['enough_space']:
        raise OSError(
            errno.ENOSPC,
            f'Insufficient disk space before {context}: need at least '
            f'{status["required_free_bytes"]:,} free bytes (download {int(download_bytes):,}, '
            f'safety factor {status["disk_safety_factor"]}, reserve '
            f'{status["disk_reserve_bytes"]:,}); only {status["free_bytes"]:,} bytes are free.',
            status['cache_root'],
        )
    return status


def get_path(cache, directory, file_name):
    info = _asset_file_info(cache, directory, file_name)
    # A local/offline cache cannot download anything. Let the local cache/path
    # check below produce a precise missing-file error instead of misleadingly
    # asking for download consent or an unknown remote size.
    if not info.get('cached', False) and not getattr(cache, '_local', False):
        _assert_download_space(
            cache,
            info.get('download_required_bytes'),
            f'downloading {directory}/{file_name}',
            unknown_files=[f'{directory}/{file_name}'],
        )
    policy = getattr(cache, '_importabcatlas_network_policy', {
        'network_retries': 0, 'network_backoff_seconds': 0,
    })
    path = _retry_call(
        lambda: cache.get_file_path(directory=directory, file_name=file_name),
        policy,
        f'downloading {directory}/{file_name}')
    if isinstance(path, dict):
        path = path['local_path']
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'File unavailable in this cache: {directory}/{file_name}')
    return path


def read_csv(path, wanted=None, preview_n=None):
    """Preserve huge numeric-looking identifiers and optionally stream by key."""
    # MapMyCells CSVs can start with comment/provenance lines. Skip only a
    # leading block; a global comment='#' would corrupt unquoted hex colors.
    leading=0
    with open(path,encoding='utf-8-sig') as stream:
        for line in stream:
            if line.startswith('#'):
                leading+=1
            else:
                break
    header = pd.read_csv(path,nrows=0,skiprows=leading)
    types = {c:'string' for c in header.columns if KEY_PATTERN.search(c)}
    kw = dict(dtype=types,low_memory=False,keep_default_na=False,na_values=[''],skiprows=leading)
    if preview_n is not None:
        positive_int(preview_n,'preview_n')
    if not wanted:
        return pd.read_csv(path,nrows=preview_n,**kw)
    missing = set(wanted)-set(header.columns)
    if missing:
        raise ValueError(f'Missing selection columns in {path}: {missing}')
    wanted = {k:set(map(str,as_list(v))) for k,v in wanted.items()}
    limit = positive_int(preview_n, 'preview_n') if preview_n is not None else None
    kept=[]; collected=0
    with pd.read_csv(path,chunksize=100000,**kw) as reader:
        for chunk in reader:
            pred=np.ones(len(chunk),dtype=bool)
            for key, vals in wanted.items():
                pred &= chunk[key].astype('string').isin(vals).to_numpy(dtype=bool,na_value=False)
            matched=chunk.loc[pred]
            if limit is not None:
                matched=matched.head(max(0,limit-collected))
            if len(matched):
                kept.append(matched); collected += len(matched)
            if limit is not None and collected >= limit:
                break
    return pd.concat(kept,ignore_index=True) if kept else header


def table(cache,directory,file_name,wanted=None,preview_n=None):
    return read_csv(get_path(cache,directory,file_name),wanted,preview_n)


def unique_gene_identifiers(values):
    """Create stable unique feature IDs without aggregating duplicate rows."""
    original = [str(value) for value in values]
    counts = pd.Series(original, dtype='string').value_counts(dropna=False).to_dict()
    used = set(original)
    seen = {}
    output = []
    for value in original:
        if counts.get(value, 0) <= 1:
            output.append(value)
            continue
        seen[value] = seen.get(value, 0) + 1
        index = seen[value]
        candidate = f'{value}.iabdup{index}'
        while candidate in used:
            index += 1
            candidate = f'{value}.iabdup{index}'
        used.add(candidate)
        output.append(candidate)
    duplicate_groups = {key: int(value) for key, value in counts.items() if value > 1}
    return output, duplicate_groups


def normalize_gene_metadata(frame, context):
    if 'gene_identifier' not in frame:
        raise ValueError(f"{context}: required column 'gene_identifier' missing.")
    if frame['gene_identifier'].isna().any() or frame['gene_identifier'].astype(str).eq('').any():
        raise ValueError(f'{context}: missing gene_identifier.')
    out = frame.copy()
    original = out['gene_identifier'].astype(str).tolist()
    unique, duplicate_groups = unique_gene_identifiers(original)
    if duplicate_groups:
        if 'source_gene_identifier' in out.columns:
            raise ValueError(f'{context}: source_gene_identifier already exists; cannot safely disambiguate duplicates.')
        original_series = pd.Series(original, dtype='string')
        # Preserve both the original identifier and its deterministic occurrence
        # within the duplicate group.  Some native HMBA H5AD files expose the
        # same duplicated stable identifier in ``var/gene_identifier`` while
        # using a different unique dataframe index.  The occurrence metadata
        # lets the package align those rows one-to-one without aggregation.
        occurrence = original_series.groupby(original_series, sort=False).cumcount() + 1
        group_size = original_series.map(original_series.value_counts(dropna=False))
        insert_at = out.columns.get_loc('gene_identifier') + 1
        out.insert(insert_at, 'source_gene_identifier', original_series)
        out.insert(insert_at + 1, 'source_gene_occurrence', occurrence.astype('Int64'))
        out.insert(insert_at + 2, 'source_gene_group_size', group_size.astype('Int64'))
        out['gene_identifier'] = pd.Series(unique, dtype='string')
    else:
        out['gene_identifier'] = out['gene_identifier'].astype('string')
    validate_ids(out, 'gene_identifier', context)
    return out, {
        'applied': bool(duplicate_groups),
        'suffix': '.iabdupN' if duplicate_groups else None,
        'duplicate_groups': duplicate_groups,
        'n_duplicate_source_identifiers': len(duplicate_groups),
        'n_disambiguated_rows': int(sum(duplicate_groups.values())),
        'occurrence_alignment': bool(duplicate_groups),
    }


def validate_ids(frame, key, context):
    if key not in frame:
        raise ValueError(f'{context}: required column {key!r} missing.')
    if frame[key].isna().any() or frame[key].astype(str).eq('').any():
        raise ValueError(f'{context}: missing {key}.')
    frame[key] = frame[key].astype('string')
    if frame[key].duplicated().any():
        examples=frame.loc[frame[key].duplicated(),key].head(3).tolist()
        raise ValueError(f'{context}: duplicate {key}: {examples}; no silent aggregation allowed.')
    return frame


def safe_join(left,right,key,suffix):
    """A left, many-to-one join; retain order, cell count and collisions."""
    if key not in left or key not in right:
        raise ValueError(f'{suffix}: join key {key} missing.')
    left=left.copy();right=right.copy()
    left[key]=left[key].astype('string');right[key]=right[key].astype('string')
    right=right.drop_duplicates()
    if right[key].isna().any() or right[key].duplicated().any():
        raise ValueError(f'{suffix}: right-hand {key} is not a non-missing unique key; join would multiply cells.')
    overlap=(set(left)&set(right))-{key}
    rename={c:f'{c}__{suffix}' for c in overlap}
    if set(rename.values()) & (set(left)|set(right)):
        raise ValueError(f'{suffix}: ambiguous annotation column suffix.')
    right=right.rename(columns=rename)
    out=left.merge(right,on=key,how='left',sort=False,validate='many_to_one')
    if len(out)!=len(left):
        raise AssertionError('Join changed the number of cells.')
    return out


def taxonomy_pivot(cache,directory):
    names=files(cache,directory)
    pivot='cluster_to_cluster_annotation_membership_pivoted'
    if pivot in names:
        result=table(cache,directory,pivot)
        return validate_ids(result,'cluster_alias',directory)
    long_name='cluster_to_cluster_annotation_membership'
    if long_name not in names:
        raise ValueError(f'{directory}: taxonomy membership file missing.')
    m=table(cache,directory,long_name)
    if 'cluster_annotation_term_name' not in m or 'cluster_annotation_term_set_name' not in m:
        if 'cluster_annotation_term' not in names:
            raise ValueError(f'{directory}: term names absent and no term table available.')
        term=table(cache,directory,'cluster_annotation_term')
        term=term.rename(columns={'label':'cluster_annotation_term_label','name':'cluster_annotation_term_name'})
        m=safe_join(m,term,'cluster_annotation_term_label','term')
        if 'cluster_annotation_term_set_name' not in m:
            sets=table(cache,directory,'cluster_annotation_term_set')
            sets=sets.rename(columns={'label':'cluster_annotation_term_set_label','name':'cluster_annotation_term_set_name'})
            m=safe_join(m,sets,'cluster_annotation_term_set_label','termset')
    keys=['cluster_alias','cluster_annotation_term_set_name']
    for key in keys+['cluster_annotation_term_name']:
        if key not in m:
            raise ValueError(f'{directory}: missing taxonomy field {key}.')
    m['cluster_alias']=m['cluster_alias'].astype('string')
    if m[keys].isna().any().any():
        raise ValueError(f'{directory}: missing taxonomy keys.')
    conflicts=m.groupby(keys)['cluster_annotation_term_name'].nunique(dropna=False)
    if (conflicts>1).any():
        raise ValueError(f'{directory}: conflicting cluster annotations.')
    result=m.pivot_table(index='cluster_alias',columns=keys[1],values='cluster_annotation_term_name',aggfunc='first')
    result.columns.name=None
    return result.reset_index()


def add_annotations(cache,spec,cell,used):
    d=spec['directory']; names=files(cache,d)
    # Preserve all donor/disease variables; no clinical interpretation is inferred.
    for fn,key in [('library','library_label'),('sample','sample_label'),('donor','donor_label'),
                   ('disease','donor_label'),('specimen_metadata','brain_section_label')]:
        if fn in names and key in cell:
            t=table(cache,d,fn,wanted={key:cell[key].dropna().astype(str).unique().tolist()})
            cell=safe_join(cell,t,key,fn);used.append([d,fn])
    if spec.get('special')=='aging':
        # Native Aging annotations are authoritative for cluster identity/age bias.
        # The separate WMB-taxonomy directory provides reference mappings.
        ad=spec.get('mapping_directory') or exact_directory('Zeng-Aging-Mouse-WMB-taxonomy',cache.list_directories)
        for td,fn,key in [(d,'cell_cluster_annotations','cell_label'),
                          (ad,'cell_cross_mapping_annotations','cell_label')]:
            if fn in files(cache,td) and key in cell:
                t=table(cache,td,fn,wanted={key:cell[key].dropna().astype(str).unique().tolist()})
                cell=safe_join(cell,t,key,fn);used.append([td,fn])
        return cell
    if spec.get('special')=='pmdbs':
        # Native PMDBS cell metadata and its own human gene table are primary.
        # SEA-AD/WHB MapMyCells tables are optional large reference mappings.
        return cell
    tax=spec.get('taxonomy')
    if tax:
        membership_dir=next((td for td in [d,tax] if 'cell_to_cluster_membership' in files(cache,td)),None)
        if membership_dir:
            t=table(cache,membership_dir,'cell_to_cluster_membership',wanted={'cell_label':cell.cell_label.tolist()})
            # Membership from the selected taxonomy is authoritative, but retain
            # any original alias under a distinct provenance-bearing name.
            if 'cluster_alias' in cell:
                cell=cell.rename(columns={'cluster_alias':'cluster_alias_original'})
            cell=safe_join(cell,t,'cell_label','membership');used.append([membership_dir,'cell_to_cluster_membership'])
        if 'cluster_alias' not in cell:
            raise ValueError(f'{d}: no cluster_alias or applicable membership table. '
                             'Use annotations=FALSE to inspect metadata; do not invent assignments.')
        t=taxonomy_pivot(cache,tax)
        cell=safe_join(cell,t,'cluster_alias','taxonomy');used.append([tax,'taxonomy annotations'])
    if spec.get('spatial'):
        for fn in ('slab_plane_coordinates','cell_anatomical_annotations'):
            if fn in names:
                t=table(cache,d,fn,wanted={'cell_label':cell.cell_label.tolist()})
                cell=safe_join(cell,t,'cell_label',fn);used.append([d,fn])
    return cell


def filter_cells(cell,filters=None,cell_ids=None,n_cells=None,seed=1):
    out=cell.copy()
    for column, values in (filters or {}).items():
        if column not in out:
            raise ValueError(f'Filter column {column!r} absent. Inspect names(cell_metadata).')
        vals=[str(v) for v in as_list(values)]
        out=out.loc[out[column].astype('string').isin(vals)]
    if cell_ids is not None:
        ids=[str(x) for x in as_list(cell_ids)]
        if len(set(ids)) != len(ids):
            raise ValueError('cell_ids contains duplicates.')
        indexed=out.set_index('cell_label',drop=False)
        missing=set(ids)-set(indexed.index)
        if missing:
            raise ValueError(f'{len(missing)} requested cells absent after filtering: {sorted(missing)[:3]}')
        out=indexed.loc[ids].reset_index(drop=True)
    if out.empty:
        raise ValueError('No cells left after filtering; inspect values and case.')
    if n_cells is not None:
        positive_int(n_cells,'n_cells')
        if n_cells>len(out):
            raise ValueError(f'Requested {n_cells} cells, but only {len(out)} qualify.')
        # Sort first for stable input-independent sampling; record IDs in provenance.
        out=out.sort_values('cell_label',kind='stable').sample(n=n_cells,random_state=int(seed))
    return out.reset_index(drop=True)


def load_data(cache,dataset,filters=None,cell_ids=None,n_cells=None,seed=1,
              preview_n=None,annotations=True,overrides=None):
    if preview_n is not None and cell_ids is not None:
        raise ValueError('preview_n and cell_ids cannot be combined.')
    spec=resolve_spec(cache,dataset,overrides)
    used=[[spec['directory'],spec['cell_file']]]
    wanted={'cell_label':as_list(cell_ids)} if cell_ids is not None else None
    if spec.get('species_filter'):
        species = spec['species_filter']
        # Selecting a preview before species filtering is incorrect whenever the
        # first rows belong to another species. Resolve species donors first,
        # then stream only their cells from the shared cell table.
        donor_dir = spec['directory']
        donor_files = files(cache, donor_dir)
        if 'donor' not in donor_files:
            aligned = next((d for d in cache.list_directories
                            if d.casefold() == 'hmba-10xmultiome-bg-aligned'), None)
            if aligned and 'donor' in files(cache, aligned):
                donor_dir = aligned
            else:
                raise ValueError(f'{dataset}: no donor metadata for species filtering.')
        donors = table(cache, donor_dir, 'donor')
        required = {'donor_label', 'species_genus'}
        if not required.issubset(donors.columns):
            raise ValueError(f'{dataset}: donor metadata lacks donor_label/species_genus.')
        selected_donors = donors.loc[
            donors.species_genus.astype(str).str.casefold() == species.casefold(),
            'donor_label'
        ].astype(str).tolist()
        if not selected_donors:
            available_species = sorted(donors.species_genus.dropna().astype(str).unique().tolist())
            raise ValueError(f'{dataset}: no {species} donors in donor metadata; '
                             f'available species_genus values: {available_species}.')
        species_wanted = {'donor_label': selected_donors}
        if wanted is not None:
            species_wanted['cell_label'] = wanted['cell_label']
        cell = table(cache, spec['directory'], spec['cell_file'],
                     wanted=species_wanted, preview_n=preview_n)
        if cell.empty:
            raise ValueError(f'{dataset}: no {species} cells found after species filtering.')
        used.append([donor_dir, 'donor'])
    else:
        cell=table(cache,spec['directory'],spec['cell_file'],wanted=wanted,preview_n=preview_n)
    cell=validate_ids(cell,'cell_label',dataset)
    if annotations:
        cell=add_annotations(cache,spec,cell,used)
    cell=filter_cells(cell,filters,cell_ids,n_cells,seed)
    gene=table(cache,spec['gene_directory'],spec['gene_file'])
    gene,gene_id_normalization=normalize_gene_metadata(gene,spec['gene_directory'])
    if 'gene_symbol' not in gene:
        gene['gene_symbol']=gene.gene_identifier
    used.append([spec['gene_directory'],spec['gene_file']])
    provenance={'package_version':VERSION,'dataset':dataset,'spec':spec,
                'manifest':str(cache.current_manifest),'metadata_files':used,
                'preview_only':preview_n is not None,'n_cells':len(cell),
                'selection_seed':int(seed),'filters':filters or {},
                'cell_id_sha256':hashlib.sha256('\n'.join(cell.cell_label).encode()).hexdigest(),
                'gene_identifier_normalization':gene_id_normalization}
    if 'gene_universe_note' in spec:
        provenance['gene_universe_note']=spec['gene_universe_note']
    return cell,gene,provenance


def resolve_genes(gene,genes=None):
    validate_ids(gene,'gene_identifier','gene metadata')
    if genes is None:
        return gene.copy().reset_index(drop=True)
    query=list(map(str,as_list(genes)))
    if not query or len(set(query))!=len(query):
        raise ValueError('genes must be a nonempty vector without duplicate requests.')
    ids=gene.gene_identifier.astype(str)
    source=(gene.source_gene_identifier.astype(str)
            if 'source_gene_identifier' in gene else None)
    symbols=(gene.gene_symbol.astype(str)
             if 'gene_symbol' in gene else pd.Series([],dtype='string'))
    selected=[]
    for q in query:
        matches=gene.index[ids.eq(q)].tolist()
        matched_kind='gene_identifier'
        if not matches and source is not None:
            matches=gene.index[source.eq(q)].tolist()
            matched_kind='source_gene_identifier'
        if not matches and len(symbols):
            matches=gene.index[symbols.eq(q)].tolist()
            matched_kind='gene_symbol'
        if len(matches)!=1:
            choices=[]
            if matches:
                choices=gene.loc[matches,'gene_identifier'].astype(str).head(8).tolist()
            detail=(f' Available disambiguated identifiers: {choices}.' if choices else '')
            raise ValueError(f'Gene {q!r} matches {len(matches)} identifiers/rows via {matched_kind}; '
                             f'use a unique gene_identifier.{detail} Absent genes are not replaced by zeros.')
        selected.append(matches[0])
    out=gene.loc[selected].reset_index(drop=True)
    if out.gene_identifier.duplicated().any():
        raise ValueError('Multiple requests resolve to the same gene identifier.')
    return out


def resolve_data_type(spec, data_type=None):
    """Resolve and validate an expression representation for one route.

    ``available_data_types`` is curated from the pinned Allen manifest.  A
    representation that is not published upstream is rejected before any
    expression object is downloaded.
    """
    available = [str(x) for x in as_list(spec.get('available_data_types', ['raw', 'log2']))]
    if not available:
        raise ValueError('This dataset route declares no accessible expression representation.')
    if data_type is None or str(data_type).strip() == '':
        data_type = spec.get('default_data_type', available[0])
    data_type = str(data_type)
    if data_type not in ('raw', 'log2'):
        raise ValueError("data_type must be one of: 'raw', 'log2'.")
    if spec.get('imputed') and data_type != 'log2':
        raise ValueError(
            'Imputed MERFISH has predicted log2 expression, not measured raw counts. '
            'Use data_type="log2".'
        )
    if data_type not in available:
        directory = spec.get('directory', 'This dataset')
        raise ValueError(
            f'{directory} does not publish data_type={data_type!r}. '
            f'Available data_type values: {available}. '
            f'The route default is {spec.get("default_data_type", available[0])!r}.'
        )
    return data_type


def expression_plan(cache,spec,cell,data_type=None,file_map=None):
    data_type = resolve_data_type(spec, data_type)
    known=list(cache.list_directories)
    preferred=[exact_directory(d,known) for d in spec['expression_directories']]
    all_inventory=[]
    for d in preferred:
        all_inventory.extend((d,name) for name in files(cache,d,'expression_matrices'))
    inventory=[item for item in all_inventory
               if item[1].casefold().endswith('/'+data_type.casefold())]
    if not inventory:
        available=sorted({name.rsplit('/',1)[-1] for _,name in all_inventory if '/' in name})
        default=spec.get('default_data_type','raw')
        raise ValueError(
            f'No {data_type} expression file is available in the configured directories '
            f'{preferred}. Available data_type values: {available or "none"}. '
            f'This route defaults to data_type={default!r}.')
    labels = cell['feature_matrix_label'].astype('string') if 'feature_matrix_label' in cell else pd.Series([spec['directory']]*len(cell))
    if labels.isna().any():
        raise ValueError('Missing feature_matrix_label in selected cells.')
    plans=[]
    for label in labels.unique():
        ids=cell.loc[labels.to_numpy()==label,'cell_label'].tolist()
        if file_map and label in file_map:
            item=file_map[label]
            mapped_directory=exact_directory(item['directory'],known)
            pair=(mapped_directory,item['file_name'])
            mapped_files=files(cache,mapped_directory,'expression_matrices')
            if item['file_name'] not in mapped_files:
                raise ValueError(f'Explicit file_map entry is not in the manifest: {pair}')
            if not item['file_name'].casefold().endswith('/'+data_type.casefold()):
                raise ValueError(f'Explicit file_map entry {pair} is not data_type={data_type!r}.')
        else:
            exact=[x for x in inventory if x[1].casefold()==f'{label}/{data_type}'.casefold()]
            matches=[x for x in exact if x[0] in preferred]
            if spec.get('imputed'):
                matches=[x for x in inventory if x[0] in preferred]
            if not matches:
                matches=[x for x in inventory if x[0] in preferred]
            if len(matches)!=1:
                raise ValueError(f'Cannot uniquely route feature_matrix_label={label!r} for '
                                 f'data_type={data_type!r}: {matches[:8]}. '
                                 f'Available expression files: {[name for _,name in all_inventory][:16]}. '
                                 'Inspect list_files(); supply file_map only after verifying the manifest.')
            pair=matches[0]
        plans.append({'directory':pair[0],'file_name':pair[1],'feature_matrix_label':str(label),
                      'n_cells':len(ids),'cell_ids':ids})
    return plans




def _annotate_expression_plans(cache, plans):
    info_by_file = {}
    annotated = []
    for item in plans:
        key = (item['directory'], item['file_name'])
        if key not in info_by_file:
            info_by_file[key] = _asset_file_info(cache, *key)
        annotated.append(dict(item) | dict(info_by_file[key]))
    return annotated


def make_fetch_plan(cache, spec, cell, gene, data_type=None, file_map=None):
    annotated = _annotate_expression_plans(
        cache, expression_plan(cache, spec, cell, data_type, file_map)
    )
    public = [{key: value for key, value in item.items() if key != 'cell_ids'}
              for item in annotated]
    # Count each source file once even if several cell groups point to it.
    unique = {(f['directory'], f['file_name']): f for f in public}
    known = all(f['file_bytes'] is not None for f in unique.values())
    total = sum(f['file_bytes'] for f in unique.values()) if known else None
    required_known = all(f['download_required_bytes'] is not None for f in unique.values())
    required = sum(f['download_required_bytes'] for f in unique.values()) if required_known else None
    return public, total, required

def h5_strings(dataset):
    """Decode H5AD strings or categorical arrays without numeric coercion."""
    if isinstance(dataset,h5py.Group):
        if 'categories' in dataset and 'codes' in dataset:
            categories=h5_strings(dataset['categories']);codes=np.asarray(dataset['codes'])
            return np.asarray([categories[i] if i>=0 else '' for i in codes],dtype=object)
        raise ValueError('Unsupported H5AD string encoding.')
    values=np.asarray(dataset)
    return np.asarray([x.decode('utf-8') if isinstance(x,(bytes,np.bytes_)) else str(x) for x in values],dtype=object)


def index_node(group):
    key=group.attrs.get('_index','_index')
    if isinstance(key,bytes):
        key=key.decode()
    if key not in group:
        if 'index' in group:
            key='index'
        else:
            raise ValueError('H5AD dataframe index is missing.')
    return group[key]


def locate_obs(obs,cell,override=None):
    names=h5_strings(index_node(obs))
    if len(set(names))!=len(names):
        raise ValueError('H5AD obs names are not unique.')
    lookup=pd.Index(names)
    columns=[override] if override else ['cell_label','exp_component_name']
    for col in columns:
        if col and col in cell and not cell[col].isna().any():
            wanted=cell[col].astype(str).tolist()
            if len(set(wanted))!=len(wanted):
                continue
            pos=lookup.get_indexer(wanted)
            if np.all(pos>=0):
                return pos,col
    # Some files store cell_label as a column rather than the dataframe index.
    if not override and 'cell_label' in obs:
        labels=h5_strings(obs['cell_label'])
        if len(set(labels))==len(labels):
            pos=pd.Index(labels).get_indexer(cell.cell_label.astype(str))
            if np.all(pos>=0):
                return pos,'obs/cell_label'
    raise ValueError('Not all selected cells can be matched exactly in H5AD. '
                     'cell_label and exp_component_name were checked; use obs_column for an explicit key. '
                     'Barcode suffixes are never stripped.')


STABLE_GENE_ALIAS_FIELDS = (
    'gene_identifier',
    'source_gene_identifier',
    'ensembl_gene_identifier',
    'ncbi_gene_identifier',
    'mapped_ncbi_identifier',
    'transcript_identifier',
)

SYMBOL_GENE_ALIAS_FIELDS = (
    'gene_symbol',
    'ncbi_gene_symbol',
)

GENE_ALIAS_FIELDS = STABLE_GENE_ALIAS_FIELDS + SYMBOL_GENE_ALIAS_FIELDS


def _clean_text(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    return text if text and text.casefold() not in {'nan', 'none', '<na>'} else None


def _identifier_variants(value):
    """Return conservative equivalent spellings for stable identifiers."""
    text = _clean_text(value)
    if text is None:
        return ()
    variants = [text]
    upper = text.upper()
    for prefix in ('ENSEMBL:', 'NCBIGENE:'):
        if upper.startswith(prefix):
            stripped = text[len(prefix):]
            if stripped:
                variants.append(stripped)
    # Preserve order while removing duplicates.
    return tuple(dict.fromkeys(variants))


def _lookup_positions(values, identifier_variants=False):
    """Return a string -> tuple(position) lookup without numeric coercion."""
    out = {}
    for index, value in enumerate(values):
        keys = _identifier_variants(value) if identifier_variants else (_clean_text(value),)
        for key in keys:
            if key is not None:
                out.setdefault(key, []).append(index)
    return {
        key: tuple(dict.fromkeys(value))
        for key, value in out.items()
    }


def _h5ad_var_aliases(var):
    """Read the H5AD feature axis and lightweight identifier columns.

    AnnData files are not fully uniform across ABC Atlas modalities: some use a
    stable identifier as ``var_names`` while others retain useful identifiers or
    symbols in ``var`` columns.  We read only identifier-like columns and never
    materialize expression here.
    """
    raw_index = h5_strings(index_node(var))
    unique_index, duplicate_groups = unique_gene_identifiers(raw_index)
    columns = {}
    for field in GENE_ALIAS_FIELDS:
        if field in var:
            try:
                values = h5_strings(var[field])
            except (TypeError, ValueError, KeyError):
                continue
            if len(values) == len(raw_index):
                columns[field] = values
    return {
        'raw_index': np.asarray(raw_index, dtype=object),
        'unique_index': np.asarray(unique_index, dtype=object),
        'columns': columns,
        'duplicate_groups': duplicate_groups,
    }


def _integer_or_none(value):
    """Return a positive integer from a pandas/scalar value, otherwise None."""
    if value is None or value is pd.NA:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    try:
        integer = int(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return integer if integer > 0 else None


def _duplicate_occurrence_from_row(row, values):
    """Return (occurrence, group_size) for a disambiguated source identifier."""
    occurrence = _integer_or_none(row.get('source_gene_occurrence'))
    group_size = _integer_or_none(row.get('source_gene_group_size'))
    if occurrence is None:
        identifier = values.get('gene_identifier')
        if identifier:
            match = re.search(r'\.iabdup([1-9][0-9]*)$', identifier)
            if match:
                occurrence = int(match.group(1))
    return occurrence, group_size


def _feature_method_priority(method):
    """Rank feature-matching evidence from stable identifiers to symbols."""
    if method is None:
        return 999
    label = str(method)
    if label == 'var_names/unique':
        return 0
    if label.startswith('var_names/source'):
        return 1
    if 'gene_identifier' in label or 'source_gene_identifier' in label:
        return 2
    if 'mapped_ncbi_identifier' in label:
        return 3
    if 'transcript_identifier' in label:
        return 4
    if 'gene_symbol' in label:
        return 20
    return 50


def _match_gene_rows_to_h5ad(var, gene, require_all=True):
    """Map gene metadata rows to exact H5AD feature positions.

    Matching is deterministic and tiered.  Stable identifiers are preferred to
    symbols.  An ambiguous alias does not immediately fail: later, stronger
    aliases are still considered.  For source identifiers duplicated both in
    metadata and in the H5AD, the package uses the preserved occurrence number
    only when the duplicate group sizes agree exactly.  No rows are aggregated.
    """
    gene = gene.copy().reset_index(drop=True)
    validate_ids(gene, 'gene_identifier', 'gene metadata')
    axis = _h5ad_var_aliases(var)
    unique_lookup = _lookup_positions(axis['unique_index'])
    raw_lookup = _lookup_positions(axis['raw_index'])
    raw_identifier_lookup = _lookup_positions(
        axis['raw_index'], identifier_variants=True
    )
    column_lookups = {
        field: _lookup_positions(
            values,
            identifier_variants=field in STABLE_GENE_ALIAS_FIELDS,
        )
        for field, values in axis['columns'].items()
    }

    positions = np.full(len(gene), -1, dtype=np.int64)
    matched_by = [None] * len(gene)
    unresolved_ambiguities = []

    def lookup_hit(lookup, value):
        if value is None:
            return None
        hits = []
        for variant in _identifier_variants(value):
            hits.extend(lookup.get(variant, ()))
        hits = tuple(dict.fromkeys(int(x) for x in hits))
        if len(hits) == 1:
            return int(hits[0])
        if len(hits) > 1:
            return hits
        return None

    for row_index, row in gene.iterrows():
        values = {
            field: _clean_text(row[field]) if field in gene.columns else None
            for field in GENE_ALIAS_FIELDS
        }
        tiers = [
            ('var_names/unique', unique_lookup, values['gene_identifier']),
            ('var_names/source', raw_identifier_lookup, values['source_gene_identifier']),
        ]

        # Stable metadata identifiers may be stored either as var_names or in a
        # differently named AnnData var column.  Native HMBA Macaque metadata,
        # for example, contains NCBI and Ensembl identifiers simultaneously.
        # Try every stable identifier before considering symbols.
        for metadata_field in STABLE_GENE_ALIAS_FIELDS:
            value = values.get(metadata_field)
            if value is None:
                continue
            if metadata_field not in {'gene_identifier', 'source_gene_identifier'}:
                tiers.append((
                    f'var_names/{metadata_field}',
                    raw_identifier_lookup,
                    value,
                ))
            for h5ad_field in STABLE_GENE_ALIAS_FIELDS:
                if h5ad_field not in column_lookups:
                    continue
                label = (
                    f'var/{h5ad_field}'
                    if h5ad_field == metadata_field
                    else f'var/{h5ad_field}<-{metadata_field}'
                )
                tiers.append((label, column_lookups[h5ad_field], value))

        # Symbols are deliberately last because they are not guaranteed unique.
        for metadata_field in SYMBOL_GENE_ALIAS_FIELDS:
            value = values.get(metadata_field)
            if value is None:
                continue
            tiers.append((
                f'var_names/{metadata_field}',
                raw_lookup,
                value,
            ))
            for h5ad_field in SYMBOL_GENE_ALIAS_FIELDS:
                if h5ad_field not in column_lookups:
                    continue
                label = (
                    f'var/{h5ad_field}'
                    if h5ad_field == metadata_field
                    else f'var/{h5ad_field}<-{metadata_field}'
                )
                tiers.append((label, column_lookups[h5ad_field], value))

        ambiguous_hits = []
        for label, lookup, value in tiers:
            hit = lookup_hit(lookup, value)
            if isinstance(hit, tuple):
                ambiguous_hits.append((label, value, hit))
                # Keep looking: a transcript or another stable column can still
                # resolve a symbol or source-ID ambiguity uniquely.
                continue
            if hit is not None:
                positions[row_index] = hit
                matched_by[row_index] = label
                break

        if positions[row_index] >= 0:
            continue

        # If multiple aliases are ambiguous, their intersection can still be a
        # unique feature.  This is stronger than occurrence-based fallback.
        if ambiguous_hits:
            candidate_sets = [set(hits) for _, _, hits in ambiguous_hits]
            intersection = set.intersection(*candidate_sets) if candidate_sets else set()
            if len(intersection) == 1:
                positions[row_index] = int(next(iter(intersection)))
                matched_by[row_index] = 'ambiguous_alias_intersection'
                continue

        occurrence, group_size = _duplicate_occurrence_from_row(row, values)
        if occurrence is not None:
            # Occurrence alignment is permitted only for a duplicated stable
            # source identifier, never for a symbol-only ambiguity.
            candidates = sorted(
                (
                    item for item in ambiguous_hits
                    if any(field in item[0] for field in STABLE_GENE_ALIAS_FIELDS)
                    or item[0] == 'var_names/source'
                ),
                key=lambda item: _feature_method_priority(item[0]),
            )
            for label, value, hits in candidates:
                if group_size is not None and len(hits) != group_size:
                    continue
                if occurrence <= len(hits):
                    positions[row_index] = int(hits[occurrence - 1])
                    matched_by[row_index] = f'{label}/occurrence'
                    break

        if positions[row_index] < 0 and ambiguous_hits:
            best = sorted(
                ambiguous_hits,
                key=lambda item: (_feature_method_priority(item[0]), item[0]),
            )[0]
            unresolved_ambiguities.append({
                'gene_identifier': values['gene_identifier'],
                'matched_by': best[0],
                'value': best[1],
                'positions': list(best[2][:8]),
            })

    if unresolved_ambiguities and require_all:
        first = unresolved_ambiguities[0]
        raise ValueError(
            f"Gene {first['gene_identifier']!r} is ambiguous in the H5AD via "
            f"{first['matched_by']}={first['value']!r}; matching positions "
            f"include {first['positions']}. Use a unique source identifier or "
            'request genes=NULL so the authoritative H5AD axis can drop unresolved '
            'metadata aliases and retain one row per measured feature.'
        )

    missing_mask = positions < 0
    if require_all and np.any(missing_mask):
        missing = gene.loc[missing_mask, 'gene_identifier'].astype(str).tolist()
        raise ValueError(
            f'{len(missing)} selected genes are unmeasured/absent in this H5AD: '
            f'{missing[:5]}. Gene metadata can be a superset of the measured '
            'feature axis; use genes=NULL to request the authoritative measured '
            'gene intersection, or choose identifiers returned by a successful fetch.'
        )
    return positions, matched_by, axis

def _canonical_metadata_rows_for_axis(gene, positions, matched_by):
    """Choose at most one metadata row for each H5AD feature position.

    Metadata tables can contain aliases or supersets in which two rows resolve
    to one measured feature.  The H5AD feature axis is authoritative.  Stable
    identifier matches win over symbol matches; ties are resolved
    deterministically by the package's unique gene identifier and row order.
    """
    by_position = {}
    for row_index, position in enumerate(positions):
        if int(position) >= 0:
            by_position.setdefault(int(position), []).append(int(row_index))

    selected = []
    collapsed = []
    for position in sorted(by_position):
        candidates = by_position[position]
        ranked = sorted(
            candidates,
            key=lambda row_index: (
                _feature_method_priority(matched_by[row_index]),
                str(gene.iloc[row_index]['gene_identifier']),
                row_index,
            ),
        )
        chosen = ranked[0]
        selected.append(chosen)
        if len(ranked) > 1:
            collapsed.append({
                'h5ad_position': int(position),
                'kept_gene_identifier': str(gene.iloc[chosen]['gene_identifier']),
                'dropped_gene_identifiers': [
                    str(gene.iloc[index]['gene_identifier']) for index in ranked[1:]
                ],
                'kept_match_method': matched_by[chosen],
            })
    return selected, collapsed


def _measured_gene_intersection(paths, gene):
    """Return one metadata row per feature measured in every selected H5AD.

    This is used when ``genes=NULL``.  It prevents a metadata superset (for
    example WMB metadata used by imputed MERFISH) from causing a fetch failure.
    Each source is first reduced to a one-to-one metadata-to-feature mapping;
    the final result is the conservative intersection across all sources.
    """
    gene = gene.copy().reset_index(drop=True)
    validate_ids(gene, 'gene_identifier', 'gene metadata')
    present = np.ones(len(gene), dtype=bool)
    per_source = []
    total_collapsed = 0

    for path in paths:
        with h5py.File(path, 'r') as handle:
            if 'var' not in handle:
                raise ValueError(f'{path}: H5AD lacks var.')
            positions, matched_by, axis = _match_gene_rows_to_h5ad(
                handle['var'], gene, require_all=False
            )

        canonical_rows, collapsed = _canonical_metadata_rows_for_axis(
            gene, positions, matched_by
        )
        current = np.zeros(len(gene), dtype=bool)
        if canonical_rows:
            current[np.asarray(canonical_rows, dtype=np.int64)] = True
        present &= current
        total_collapsed += sum(len(item['dropped_gene_identifiers']) for item in collapsed)

        per_source.append({
            'path': str(path),
            'h5ad_features': int(len(axis['raw_index'])),
            'metadata_rows_matched_before_canonicalization': int((positions >= 0).sum()),
            'metadata_rows_selected_one_to_one': int(current.sum()),
            'metadata_rows_missing': int((positions < 0).sum()),
            'metadata_alias_rows_collapsed': int(
                sum(len(item['dropped_gene_identifiers']) for item in collapsed)
            ),
            'match_methods': dict(pd.Series(
                [matched_by[index] for index in canonical_rows
                 if matched_by[index] is not None],
                dtype='string',
            ).value_counts().astype(int)),
            'collapsed_examples': collapsed[:8],
        })

    selected = gene.loc[present].reset_index(drop=True)
    if selected.empty:
        raise ValueError(
            'No gene metadata rows are measured in every selected H5AD source. '
            'Inspect the route and source feature axes; no zero-filled union is created.'
        )

    # Verify that the final metadata selection remains one-to-one in every
    # source.  This protects dense h5py indexing and, more importantly, avoids
    # returning the same measured feature under two different metadata rows.
    for path in paths:
        with h5py.File(path, 'r') as handle:
            positions, _, _ = _match_gene_rows_to_h5ad(
                handle['var'], selected, require_all=True
            )
        if len(np.unique(positions)) != len(positions):
            raise ValueError(
                f'{path}: authoritative reconciliation still maps multiple gene '
                'metadata rows to the same H5AD feature. Inspect identifier aliases.'
            )

    return selected, {
        'mode': 'authoritative_h5ad_intersection',
        'candidate_metadata_rows': int(len(gene)),
        'selected_measured_rows': int(len(selected)),
        'dropped_unmeasured_or_alias_rows': int(len(gene) - len(selected)),
        # Backward-compatible name retained for callers that already inspect it.
        'dropped_unmeasured_rows': int(len(gene) - len(selected)),
        'metadata_alias_rows_collapsed': int(total_collapsed),
        'sources': per_source,
    }

def subset_h5ad(path, cell, gene_ids, obs_column=None, max_output_mb=1024,
                 return_feature_matches=False):
    """Read only requested rows/columns of dense, CSR or CSC H5AD X.

    ``gene_ids`` may be a vector of identifiers or a gene-metadata DataFrame.
    When metadata are supplied, stable IDs and symbols are reconciled against
    AnnData ``var_names`` and identifier columns without guessing or aggregation.
    """
    if isinstance(gene_ids, pd.DataFrame):
        selected_gene = gene_ids.copy().reset_index(drop=True)
    else:
        selected_gene = pd.DataFrame({
            'gene_identifier': pd.Series(list(map(str, gene_ids)), dtype='string')
        })
    validate_ids(selected_gene, 'gene_identifier', 'selected gene metadata')
    if not np.isfinite(float(max_output_mb)) or float(max_output_mb) <= 0:
        raise ValueError('max_output_mb must be positive and finite.')
    rows = len(cell)
    cols = len(selected_gene)
    if rows * cols * 8 > float(max_output_mb) * 1024**2:
        raise MemoryError(
            'Requested subset exceeds max_output_mb (dense-equivalent guard). '
            'Reduce cells/genes or raise it explicitly.'
        )
    with h5py.File(path, 'r') as f:
        for key in ('obs', 'var', 'X'):
            if key not in f:
                raise ValueError(f'{path}: H5AD lacks {key}.')
        ri, matched_by = locate_obs(f['obs'], cell, obs_column)
        ci, feature_methods, axis = _match_gene_rows_to_h5ad(
            f['var'], selected_gene, require_all=True)
        if len(np.unique(ci)) != len(ci):
            duplicated_positions = pd.Series(ci).duplicated(keep=False).to_numpy()
            duplicated_genes = selected_gene.loc[
                duplicated_positions, 'gene_identifier'
            ].astype(str).tolist()
            raise ValueError(
                'Selected gene identifiers do not map one-to-one to the H5AD '
                f'feature axis; duplicated feature mappings involve {duplicated_genes[:8]}. '
                'Use genes=NULL to reconcile the authoritative measured axis, then '
                'make explicit requests with the returned gene_identifier values.'
            )
        names = axis['unique_index']
        x = f['X']
        rr = []
        cc = []
        vv = []
        if isinstance(x, h5py.Dataset):
            if x.shape != (len(h5_strings(index_node(f['obs']))), len(names)):
                raise ValueError('H5AD X shape does not agree with obs/var.')
            # h5py requires a strictly increasing, duplicate-free fancy
            # index. np.unique also returns the inverse needed to restore the
            # caller's requested gene order.
            sorted_ci, inverse = np.unique(ci, return_inverse=True)
            for out_row, r in enumerate(ri):
                val = np.asarray(x[int(r), sorted_ci])[inverse].astype(np.float64)
                nz = np.flatnonzero(val)
                rr.extend([out_row] * len(nz))
                cc.extend(nz)
                vv.extend(val[nz])
        else:
            encoding = x.attrs.get('encoding-type', '')
            if isinstance(encoding, bytes):
                encoding = encoding.decode()
            shape = tuple(x.attrs.get('shape', ()))
            if shape and shape != (len(h5_strings(index_node(f['obs']))), len(names)):
                raise ValueError('Sparse X shape disagrees with obs/var.')
            if encoding == 'csr_matrix':
                remap = {}
                for output_column, source_column in enumerate(ci):
                    remap.setdefault(int(source_column), []).append(output_column)
                for out_row, r in enumerate(ri):
                    start, end = np.asarray(x['indptr'][int(r):int(r) + 2], dtype=np.int64)
                    for at in range(int(start), int(end), 65536):
                        idx = np.asarray(x['indices'][at:min(at + 65536, int(end))], dtype=np.int64)
                        val = np.asarray(x['data'][at:min(at + 65536, int(end))], dtype=np.float64)
                        for j, v in zip(idx, val):
                            if int(j) in remap and v != 0:
                                for output_column in remap[int(j)]:
                                    rr.append(out_row)
                                    cc.append(output_column)
                                    vv.append(v)
            elif encoding == 'csc_matrix':
                remap = {int(r): i for i, r in enumerate(ri)}
                for out_col, c in enumerate(ci):
                    start, end = np.asarray(x['indptr'][int(c):int(c) + 2], dtype=np.int64)
                    for at in range(int(start), int(end), 65536):
                        idx = np.asarray(x['indices'][at:min(at + 65536, int(end))], dtype=np.int64)
                        val = np.asarray(x['data'][at:min(at + 65536, int(end))], dtype=np.float64)
                        for j, v in zip(idx, val):
                            if int(j) in remap and v != 0:
                                rr.append(remap[int(j)])
                                cc.append(out_col)
                                vv.append(v)
            else:
                raise ValueError(f'Unsupported sparse H5AD encoding: {encoding!r}.')
        values = np.asarray(vv, dtype=np.float64)
        if not np.isfinite(values).all():
            raise ValueError('Selected expression contains NaN or infinite values.')
        matrix = sparse.coo_matrix((values, (rr, cc)), shape=(rows, cols)).tocsr()
        if return_feature_matches:
            return matrix, matched_by, {
                'methods': feature_methods,
                'h5ad_feature_count': int(len(names)),
            }
        return matrix, matched_by


def fetch_data(cache, spec, cell, gene, genes=None, data_type=None,
               allow_downloads=False, max_output_mb=1024, file_map=None,
               obs_column=None, counts_semantics=None, verify_reference=False):
    data_type = resolve_data_type(spec, data_type)
    cell = validate_ids(cell.copy(), 'cell_label', 'fetch_data metadata')
    explicit_gene_request = genes is not None
    selected = resolve_genes(gene.copy(), genes) if explicit_gene_request else None
    if not np.isfinite(float(max_output_mb)) or float(max_output_mb) <= 0:
        raise ValueError('max_output_mb must be finite and positive.')
    if selected is not None and len(cell) * len(selected) * 8 > float(max_output_mb) * 1024**2:
        raise MemoryError('Requested subset exceeds max_output_mb before downloading.')

    plans = _annotate_expression_plans(
        cache, expression_plan(cache, spec, cell, data_type, file_map)
    )
    unique = {(p['directory'], p['file_name']): p for p in plans}
    missing = [p for p in unique.values() if not p.get('cached', False)]
    if missing and getattr(cache, '_local', False):
        labels = ', '.join(f"{p['directory']}/{p['file_name']}" for p in missing)
        raise FileNotFoundError(
            'Offline cache lacks required expression files: ' + labels + '. '
            'Run once with offline=FALSE and allow_downloads=TRUE after reviewing plan_fetch().'
        )
    if missing and not allow_downloads:
        labels = ', '.join(f"{p['directory']}/{p['file_name']}" for p in missing)
        raise PermissionError(
            'Expression download consent is required for: ' + labels + '. '
            'Set allow_downloads=TRUE only after reviewing plan_fetch(). '
            'Already verified cached files can be read without download consent.'
        )
    unknown = [f"{p['directory']}/{p['file_name']}" for p in missing
               if p.get('download_required_bytes') is None]
    download_bytes = (None if unknown else
                      sum(int(p['download_required_bytes']) for p in missing))
    _assert_download_space(
        cache, download_bytes, 'fetch_data expression download', unknown_files=unknown
    )

    # Resolve/download every unique source once.  The paths are then reused for
    # feature-axis reconciliation and expression extraction.
    path_by_key = {}
    for key, plan in unique.items():
        path_by_key[key] = get_path(cache, plan['directory'], plan['file_name'])

    if selected is None:
        selected, gene_axis = _measured_gene_intersection(
            list(path_by_key.values()), gene.copy()
        )
    else:
        gene_axis = {
            'mode': 'explicit_strict_request',
            'candidate_metadata_rows': int(len(gene)),
            'selected_measured_rows': int(len(selected)),
            'dropped_unmeasured_rows': 0,
            'sources': [],
        }
    if len(cell) * len(selected) * 8 > float(max_output_mb) * 1024**2:
        raise MemoryError(
            'Requested measured subset exceeds max_output_mb after reconciling '
            'the authoritative H5AD feature axis.'
        )

    parts = []
    ordered_ids = []
    sources = []
    feature_match_summaries = []
    for p in plans:
        sub = cell.set_index('cell_label', drop=False).loc[p['cell_ids']].reset_index(drop=True)
        path = path_by_key[(p['directory'], p['file_name'])]
        m, matched_by, feature_matches = subset_h5ad(
            path, sub, selected, obs_column, max_output_mb,
            return_feature_matches=True
        )
        if verify_reference:
            verify_with_anndata(path, sub, selected, m, matched_by)
        parts.append(m)
        ordered_ids.extend(sub.cell_label.tolist())
        method_counts = dict(pd.Series(
            feature_matches['methods'], dtype='string'
        ).value_counts().astype(int))
        feature_match_summaries.append({
            'directory': p['directory'],
            'file_name': p['file_name'],
            'h5ad_feature_count': feature_matches['h5ad_feature_count'],
            'match_methods': method_counts,
        })
        sources.append({k: v for k, v in p.items() if k != 'cell_ids'} |
                       {'local_path': str(path), 'file_bytes': path.stat().st_size,
                        'matched_by': matched_by,
                        'feature_match_methods': method_counts})
    joined = sparse.vstack(parts, format='csr')
    order = pd.Index(ordered_ids).get_indexer(cell.cell_label)
    if np.any(order < 0) or len(set(ordered_ids)) != len(cell):
        raise AssertionError('Expression/cell alignment invariant failed.')
    matrix = joined[order, :].T.tocsc()
    integer_nonnegative = bool(
        np.all(matrix.data >= 0) and
        np.allclose(matrix.data, np.rint(matrix.data), rtol=0, atol=1e-8)
    )
    if counts_semantics is None:
        is_counts = (data_type == 'raw' and not spec.get('imputed', False)
                     and integer_nonnegative)
    else:
        is_counts = bool(counts_semantics)
        if is_counts and (data_type != 'raw' or spec.get('imputed')
                          or not integer_nonnegative):
            raise ValueError(
                'Cannot label transformed, imputed, negative or noninteger data as counts.'
            )
    assay_kind = ('counts' if is_counts else
                  ('imputed_log2' if spec.get('imputed') else data_type))
    provenance = {
        'package_version': VERSION,
        'manifest': str(cache.current_manifest),
        'sources': sources,
        'data_type': data_type,
        'assay_kind': assay_kind,
        'integer_nonnegative': integer_nonnegative,
        'reference_comparison': bool(verify_reference),
        'n_cells': matrix.shape[1],
        'n_genes': matrix.shape[0],
        'gene_request_mode': 'explicit' if explicit_gene_request else 'all_measured_intersection',
        'gene_axis_reconciliation': gene_axis,
        'feature_match_summaries': feature_match_summaries,
        'cell_id_sha256': hashlib.sha256('\n'.join(cell.cell_label).encode()).hexdigest(),
        'gene_id_sha256': hashlib.sha256(
            '\n'.join(selected.gene_identifier).encode()
        ).hexdigest(),
    }
    return matrix, cell, selected, provenance


def verify_with_anndata(path, cell, gene, actual, matched_by):
    """Independent source-value comparison, intended for small live smoke tests.

    ``gene`` may be the metadata DataFrame used by the backend or a vector of
    requested gene identifiers.  Normalize both forms to the canonical
    ``gene_identifier`` DataFrame before feature-axis reconciliation.
    """
    try:
        import anndata
    except ImportError as error:
        raise ImportError('anndata is required for verify_reference=True.') from error

    if isinstance(gene, pd.DataFrame):
        gene = gene.copy().reset_index(drop=True)
    elif isinstance(gene, pd.Series):
        gene = pd.DataFrame({
            'gene_identifier': gene.astype(str).tolist()
        })
    elif isinstance(gene, (list, tuple, np.ndarray, pd.Index)):
        gene = pd.DataFrame({
            'gene_identifier': [str(value) for value in gene]
        })
    else:
        raise TypeError(
            'gene must be a pandas DataFrame, Series, Index, list, tuple, or NumPy array.'
        )

    if 'gene_identifier' not in gene.columns:
        raise ValueError(
            "gene DataFrame must contain a 'gene_identifier' column for independent verification."
        )
    ad = anndata.read_h5ad(path, backed='r')
    try:
        axis = (pd.Index(ad.obs['cell_label'].astype(str))
                if matched_by == 'obs/cell_label'
                else pd.Index(ad.obs_names.astype(str)))
        query = cell['cell_label' if matched_by == 'obs/cell_label'
                     else matched_by].astype(str)
        ri = axis.get_indexer(query)
        with h5py.File(path, 'r') as handle:
            ci, _, _ = _match_gene_rows_to_h5ad(
                handle['var'], gene, require_all=True
            )
        if np.any(ri < 0) or np.any(ci < 0):
            raise AssertionError('Independent AnnData axis lookup disagrees with importer.')
        for i, row_index in enumerate(ri):
            row = ad[int(row_index), :].X
            expected = (row.toarray() if sparse.issparse(row)
                        else np.asarray(row)).reshape(-1)[ci]
            np.testing.assert_allclose(
                actual.getrow(i).toarray().reshape(-1), expected,
                rtol=0, atol=0,
                err_msg='Independent AnnData source-value comparison failed.'
            )
    finally:
        ad.file.close()

def environment():
    versions={}
    for package in ('numpy','pandas','scipy','h5py','abc_atlas_access'):
        try: versions[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: versions[package]=None
    try:
        dist=importlib.metadata.distribution('abc_atlas_access')
        direct=dist.read_text('direct_url.json')
        versions['abc_atlas_access_source']=json.loads(direct) if direct else None
    except importlib.metadata.PackageNotFoundError:
        pass
    return {'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
            'packages':versions,'package_version':VERSION}


def write_frame(frame,path):
    frame.to_csv(path,index=False)
    schema={c:('string' if KEY_PATTERN.search(c) else
               'numeric' if pd.api.types.is_numeric_dtype(frame[c]) and not pd.api.types.is_bool_dtype(frame[c])
               else 'logical' if pd.api.types.is_bool_dtype(frame[c]) else 'string') for c in frame}
    Path(str(path)+'.schema.json').write_text(json.dumps(schema))


def file_catalog(cache,directory=None):
    rows=[]
    directories=[exact_directory(directory,cache.list_directories)] if directory else cache.list_directories
    for d in directories:
        for kind in KINDS:
            for name in files(cache,d,kind):
                rows.append({'directory':d,'kind':kind,'file_name':name})
    return rows


def run_request(req):
    action=req['action'];out=Path(req['out_dir']);out.mkdir(parents=True,exist_ok=True)
    args=req.get('args',{})
    if action=='diagnose':return environment()
    if action=='registry':return REGISTRY
    if action=='demo':
        root=Path(__file__).resolve().parents[1]/'extdata'
        cell=read_csv(root/'demo_cells.csv');gene=read_csv(root/'demo_genes.csv')
        matrix,_=subset_h5ad(root/'demo.h5ad',cell,gene.gene_identifier.tolist())
        matrix=matrix.T.tocsc()
        write_frame(cell,out/'cell_metadata.csv');write_frame(gene,out/'gene_data.csv');mmwrite(out/'matrix.mtx',matrix)
        return {'assay_kind':'counts','synthetic_fixture':True,'n_cells':len(cell),'n_genes':len(gene)}
    cache=create_cache(
        args.pop('download_base','abc_download_root'),
        args.pop('release',DEFAULT_RELEASE),
        args.pop('offline',False),
        check_disk_space=args.pop('check_disk_space',True),
        disk_reserve_gb=args.pop('disk_reserve_gb',1.0),
        disk_safety_factor=args.pop('disk_safety_factor',1.05),
        allow_unknown_size=args.pop('allow_unknown_size',False),
        network_retries=args.pop('network_retries',5),
        network_backoff_seconds=args.pop('network_backoff_seconds',2.0),
        network_connect_timeout_seconds=args.pop('network_connect_timeout_seconds',30),
        network_read_timeout_seconds=args.pop('network_read_timeout_seconds',300),
        download_chunk_mb=args.pop('download_chunk_mb',8),
    )
    if action=='catalog':return {'manifest':cache.current_manifest,'files':file_catalog(cache,args.get('directory'))}
    if action=='download_file':
        path=get_path(cache,args['directory'],args['file_name'])
        return {'path':str(path),'manifest':cache.current_manifest}
    if action=='metadata':
        frame=table(cache,args['directory'],args['file_name'],preview_n=args.get('preview_n'))
        write_frame(frame,out/'table.csv');return {'manifest':cache.current_manifest}
    if action=='load':
        cell,gene,prov=load_data(cache,**args)
        write_frame(cell,out/'cell_metadata.csv');write_frame(gene,out/'gene_data.csv')
        return prov
    if action in ('fetch','plan'):
        cell=read_csv(args.pop('metadata_path'));gene=read_csv(args.pop('gene_path'))
        cell=filter_cells(cell,args.pop('filters',None),args.pop('cell_ids',None),args.pop('n_cells',None),args.pop('seed',1))
        spec=resolve_spec(cache,args.pop('dataset'),args.pop('overrides',None))
        if action=='plan':
            requested_genes=args.get('genes')
            selected=resolve_genes(gene,requested_genes)
            resolved_data_type=resolve_data_type(spec,args.get('data_type'))
            plans,total_bytes,download_bytes=make_fetch_plan(
                cache,spec,cell,gene,resolved_data_type,args.get('file_map'))
            disk=_disk_status(cache,download_bytes)
            gene_count_basis=(
                'metadata_candidate_upper_bound; authoritative measured intersection is resolved during fetch'
                if requested_genes is None else
                'explicit metadata selection; availability is validated against each H5AD during fetch'
            )
            unknown=[f"{f['directory']}/{f['file_name']}" for f in plans if f.get('file_bytes') is None]
            size_errors=[f.get('size_error') for f in plans if f.get('size_error')]
            warning=('No expression files were downloaded. Exact sizes are per selected file, obtained from '
                     'a verified local file, manifest attributes, or an S3 object-size request. '
                     + (('Size remains unknown for: '+', '.join(unknown)+'. ') if unknown else '')
                     + (('Size lookup details: '+' | '.join(size_errors)) if size_errors else ''))
            return {'manifest':cache.current_manifest,'files':plans,
                    'total_file_bytes':total_bytes,'total_download_bytes':download_bytes,
                    'data_type':resolved_data_type,
                    'n_cells':len(cell),'n_genes':len(selected),
                    'gene_count_basis':gene_count_basis,
                    'gene_count_exact_before_fetch':False,
                    'dense_equivalent_mb':len(cell)*len(selected)*8/1024**2,
                    'disk':disk,'unknown_size_files':unknown,'warning':warning}
        args['data_type']=resolve_data_type(spec,args.get('data_type'))
        matrix,cell,gene,prov=fetch_data(cache,spec,cell,gene,**args)
        # MatrixMarket plus CSV sidecars must also fit in the temporary output volume.
        sparse_bytes=int(matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes)
        temp_free=int(shutil.disk_usage(out).free)
        temp_required=max(10*1024**2,int(math.ceil(sparse_bytes*2.5)))
        if temp_free < temp_required:
            raise OSError(errno.ENOSPC,
                          f'Insufficient temporary-disk space: need about {temp_required:,} bytes; '
                          f'only {temp_free:,} bytes are free.',str(out))
        write_frame(cell,out/'cell_metadata.csv');write_frame(gene,out/'gene_data.csv');mmwrite(out/'matrix.mtx',matrix)
        prov['temporary_output_estimate_bytes']=temp_required
        return prov
    raise ValueError(f'Unknown action: {action}')


def main():
    if len(sys.argv)!=2:
        raise SystemExit('Usage: python abc_backend.py request.json')
    req=json.loads(Path(sys.argv[1]).read_text());out=Path(req['out_dir']);out.mkdir(parents=True,exist_ok=True)
    try:
        start=time.perf_counter();result=run_request(req)
        if isinstance(result,dict) and req['action']!='registry':result['backend_elapsed_seconds']=time.perf_counter()-start
        (out/'result.json').write_text(json.dumps(result,indent=2,default=str,allow_nan=False))
    except Exception as error:
        (out/'error.json').write_text(json.dumps({'type':type(error).__name__,'message':str(error),'traceback':traceback.format_exc()},indent=2))
        traceback.print_exc();return 1
    return 0


if __name__=='__main__':
    raise SystemExit(main())
