function sem = calc_sem(data)

    sem = std(data(:))/sqrt(length(data(:)));

end
