function sem = calcSEM(data)

    sem = std(data(:),'omitnan')/sqrt(sum(~isnan(data(:))));

end
