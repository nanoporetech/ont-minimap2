# ONT Minimap2

[![build-minimap2](https://github.com/nanoporetech/ont-minimap2/actions/workflows/build.yml/badge.svg)](https://github.com/nanoporetech/ont-minimap2/actions/workflows/build.yml)

Cross platform builds for [minimap2](https://github.com/lh3/minimap2/).

```
$ git clone --recurse-submodules https://github.com/nanoporetech/ont-minimap2.git
$ cmake -S . -B cmake-build
$ cmake --build cmake-build --config Release -j
$ ctest -C Release --test-dir cmake-build --output-on-failure
```

## Build Options

| CMake                    | Description                                     | Default |
|:-------------------------|:------------------------------------------------|---------|
| ONT_MM2_EXE              | Build the minimap2 executable                   | ON      |
| ONT_MM2_LIB              | Build the minimap2 library                      | OFF     |
