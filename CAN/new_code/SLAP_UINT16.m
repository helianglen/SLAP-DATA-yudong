function output = SLAP_UINT16(input, swap_endianness)

if swap_endianness
    output = typecast(uint8(flip(input)), 'uint16');
else
    output = typecast(uint8(input), 'uint16');
end