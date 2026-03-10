%%% get_framerate
% Searches through timing struct to find the shortest duration.
% The duration of a frame will be set as half the shortest duration.

function [frame_dur, frame_rate] = get_framerate(t)

struct_fieldnames = fieldnames(t);
struct_fieldnames(~contains(struct_fieldnames,{'dur','rate'})) = [];
struct_fieldnames(contains(struct_fieldnames,'monitor')) = [];

for i_field = 1:numel(struct_fieldnames)
    if i_field == 1

        min_dur = t.(struct_fieldnames{i_field});

    else
        if t.(struct_fieldnames{i_field}) < min_dur

            min_dur = t.(struct_fieldnames{i_field});

        end

    end
end

frame_dur = min_dur/2;
frame_rate = 1/frame_dur;

end