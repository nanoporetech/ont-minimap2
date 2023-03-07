# ONT Minimap2

Wrapper around minimap2 that build a shared library.

```
$ git clone --recurse-submodules https://github.com/nanoporetech/ont_minimap2.git
$ cmake -S . -B cmake-build
$ cmake --build cmake-build --config Release -j
$ ctest -C Release --test-dir cmake-build --output-on-failure
```

## Build Options

| CMake                    | Description                                     | Default |
|:-------------------------|:------------------------------------------------|---------|
| ONT_MM2_EXE              | Build the minimap2 executable                   | ON      |
