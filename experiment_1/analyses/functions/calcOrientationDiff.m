function diff = calcOrientationDiff(a1, a2)

    diff = mod((a1 - a2) + 90, 180) - 90;

end