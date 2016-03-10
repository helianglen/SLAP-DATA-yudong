function output = SLAP_UINT64(input, swap_endianness)

if swap_endianness
    output = typecast(uint8(flip(input)), 'uint64');
else
    output = typecast(uint8(input), 'uint64');
end