import importlib.metadata

try:
    __version__ = importlib.metadata.version("funannotate2_addons")
except importlib.metadata.PackageNotFoundError:
    # Fallback version for development
    __version__ = "25.5.24"
