#!/bin/bash

# Runs clang-format over source files

find lib -iname *.hpp -o -iname *.cpp -o -iname *.h | xargs clang-format -i

