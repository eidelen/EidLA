#!/bin/bash

# Runs clang-format over source files

find lib/ -iname *.hpp -o -iname *.cpp | xargs clang-format -i