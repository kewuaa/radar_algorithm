# requires

+ Python 3.8+
+ CMake 3.15+
+ C++ compiler support c++17

ensure python counld be found by cmake, and `nanobind` has been installed.
use `pip install nanobind` to install it, if use virtual environment, activate
it at first.

# build

```sh
# build python extension
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=`pwd`/src
cmake --build build
cmake --install build

# build python package wheel at dist/
pip install build
python -m build
# or use uv
uv build
```
