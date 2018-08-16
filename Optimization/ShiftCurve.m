function [shifted_curve] = ShiftCurve(curve_to_shift, shift_factor)

if nnz(shift_factor)
    shift_idx = abs(shift_factor);
    if shift_factor > 0
        shifted_curve = [zeros(shift_idx,1); curve_to_shift(1:end-shift_idx,1)];
%         shifted_curve = [curve_to_shift(end-shift_idx+1:end,1); curve_to_shift(1:end-shift_idx,1)];

    elseif shift_factor <= 0
        shifted_curve = [curve_to_shift(shift_idx:end,1); zeros(shift_idx-1,1)];
%         shifted_curve = [curve_to_shift(shift_idx:end,1); curve_to_shift(1:shift_idx-1,1)];

    end
else
    shifted_curve = curve_to_shift;
end


end