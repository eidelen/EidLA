#!/bin/bash

# Runs clang-format over source files

find lib/ -iname *.h -o -iname *.cpp -iname *.hpp | xargs clang-format -i