#!/bin/bash

touch extensions/python/pylatticeflow.pyx ; rm -rf build/extensions/python/ build/pylatticeflow.so ; ./waf build install --test --debug-full --prefix=$SYSROOT
