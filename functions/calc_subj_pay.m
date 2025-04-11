
function calc_subj_pay(subj_ID)

    total_hours = [];
    dollars_per_hour = 15;

    subj_pay = round(total_hours * dollars_per_hour);

    disp(['Payment: $' num2str(subj_pay)])

end
