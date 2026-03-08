import importlib.metadata

try:
    __version__ = importlib.metadata.version("funannotate2_addons")
except importlib.metadata.PackageNotFoundError:
    # Fallback version for development
    __version__ = "26.3.7"
