#!/bin/bash
# Print the RPC-related symbols of the embedded gfx offload image of an executable.
# usage: check_image.sh <exe>
source "$(dirname "$0")/../common.sh"
exe=$1; tmp=$(mktemp -d)
llvm-objcopy --dump-section=.llvm.offloading="$tmp/off.bin" "$exe" /dev/null || { echo "no .llvm.offloading section"; exit 1; }
# the section is an offload bundle; the device ELF starts at the first \x7fELF
off=$(grep -boa $'\x7fELF' "$tmp/off.bin" | head -1 | cut -d: -f1)
tail -c +$((off+1)) "$tmp/off.bin" > "$tmp/dev.elf"
n=$(llvm-readelf --dyn-syms "$tmp/dev.elf" 2>/dev/null | grep -c -E 'rpc|_FortranAio')
echo "$exe: $n rpc/Fortran-I/O symbols in device image"
llvm-readelf --dyn-syms "$tmp/dev.elf" 2>/dev/null | grep -E 'rpc|_FortranAio|_FortranAStop' | awk '{print "    "$NF}' | sort -u | head
rm -rf "$tmp"
