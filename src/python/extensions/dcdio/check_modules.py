import pkgutil
import sasmol

def list_modules(package):
    package_path = package.__path__
    for module_info in pkgutil.iter_modules(package_path):
        print(module_info.name)

list_modules(sasmol)
