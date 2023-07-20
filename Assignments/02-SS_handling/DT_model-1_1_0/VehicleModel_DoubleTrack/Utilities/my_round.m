function rounded_val = my_round(value,digits)

    if (value*10^digits - floor(value*10^digits) >= 0.5)
        rounded_val = ceil(value*10^digits)/10^digits;
    else
        rounded_val = floor(value*10^digits)/10^digits;
    end

end

