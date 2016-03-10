function output = SLAP_UINT32(input, swap_endianness)

if swap_endianness
    output = typecast(uint8(flip(input)), 'uint32');
else
    output = typecast(uint8(input), 'uint32');
end