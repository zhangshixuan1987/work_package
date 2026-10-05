"""Deprecated compatibility aliases; import :mod:`mov_config` instead."""

from mov_config import (  # noqa: F401
    DEFAULT_MODEL_DATA_ROOT,
    DEFAULT_MODEL_SUBDIR,
    DEFAULT_PERIOD,
    MODE_CONFIG,
    OBS_SOURCES,
    V3_RRM_MODELS as MODEL_DATASETS,
    V3_RRM_MODEL_ORDER as MODEL_ORDER,
    V3_RRM_PCMDI_MODEL_TAGS as PCMDI_MODEL_TAGS,
    v3_rrm_model_datasets as pathway_model_datasets,
)

