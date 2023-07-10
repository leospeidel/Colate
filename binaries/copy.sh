#!/bin/bash

prev_version="0.1.3"
next_version="0.1.4"

#mkdir old_versions
tar -xvf old_versions.tar
cp -r v${prev_version}_MacOSX_Intel/ old_versions
cp -r v${prev_version}_MacOSX_M1/ old_versions
cp -r v${prev_version}_x86_64_dynamic/ old_versions
cp -r v${prev_version}_x86_64_static/ old_versions

tar -rvf old_versions.tar old_versions

mv ./v${prev_version}_MacOSX_Intel/ ./v${next_version}_MacOSX_Intel/
mv ./v${prev_version}_MacOSX_M1/ ./v${next_version}_MacOSX_M1/
mv ./v${prev_version}_x86_64_dynamic/ ./v${next_version}_x86_64_dynamic/
mv ./v${prev_version}_x86_64_static/ ./v${next_version}_x86_64_static/

rm -rf old_versions
