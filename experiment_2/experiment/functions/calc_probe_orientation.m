
function probe_orientation = calc_probe_orientation(test_orientation, probe_offset)

    probe_orientation = test_orientation + probe_offset .* datasample([-1 1], length(test_orientation))';
