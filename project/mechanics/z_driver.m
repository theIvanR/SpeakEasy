function Z = z_driver(s, p)
%Z_DRIVER  Driver mechanical impedance.

    Z = p.Rms + s*p.Mms + 1 ./ (s*p.Cms);
end