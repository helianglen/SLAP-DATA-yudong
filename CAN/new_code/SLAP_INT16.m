function output = SLAP_INT16(input, swap_endianness)

if swap_endianness
    output = typecast(uint8(flip(input)), 'int16');
else
    output = typecast(uint8(input), 'int16');
end