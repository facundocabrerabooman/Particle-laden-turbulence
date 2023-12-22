function [neighbor_terminalState, neighbor_terminalState_t]  = neighborLayerTS(neighbor_layer,startTime,endTime,Nexp)

Nframemax = 1e6;
tstart = startTime(Nexp);
tend = endTime(Nexp);
tstart = tstart-Nframemax*(Nexp-1);
tend = tend-Nframemax*(Nexp-1);
tlen = tend-tstart+1;

fields = fieldnames(neighbor_layer(1,1));
neighbor_terminalState_t = struct();
for t = 1:tlen
    tt = t+tstart-1;
    for k = 1:numel(neighbor_layer(t,:))
        for i = 1:numel(fields)
            neighbor_terminalState_t(t,k).(fields{i}) = neighbor_layer(tt,k).(fields{i});
        end
    end
end


for k = 1:numel(neighbor_layer(1,:))
    for i = 1:5
        neighbor_terminalState(1,k).(fields{i}) = neighbor_terminalState_t(1,k).(fields{i});
    end
end

for k = 1:numel(neighbor_layer(1,:))
    for i = 6:9
        neighbor_terminalState(1,k).(fields{i}) = vertcat(neighbor_terminalState_t(:,k).(fields{i}));
    end
end


for k = 1:numel(neighbor_layer(1,:))
    for i = 10:25
        neighbor_terminalState(1,k).(fields{i}) = 0;
    end
end


for k = 1:numel(neighbor_layer(1,:))
    nfront = 0;
    nback  = 0;
    for j = 1:size(neighbor_terminalState_t,1)
        nnfront = numel(neighbor_terminalState_t(j,k).(fields{8}));
        nnback  = numel(neighbor_terminalState_t(j,k).(fields{9}));
        nfront  = nfront + nnfront;
        nback   = nback  + nnback;

        if nnfront ~=0
            for i = 10:13
                neighbor_terminalState(1,k).(fields{i}) = neighbor_terminalState(1,k).(fields{i}) + nnfront*neighbor_terminalState_t(j,k).(fields{i});    
            end
        end

        if nnback ~=0
            for i = 14:17
                neighbor_terminalState(1,k).(fields{i}) = neighbor_terminalState(1,k).(fields{i}) + nnback*neighbor_terminalState_t(j,k).(fields{i});
            end
        end

        if nnfront ~=0
            for i = 18:21
                neighbor_terminalState(1,k).(fields{i}) = neighbor_terminalState(1,k).(fields{i}) + nnfront*neighbor_terminalState_t(j,k).(fields{i});
            end
        end

        if nnback ~=0
            for i = 22:25
                neighbor_terminalState(1,k).(fields{i}) = neighbor_terminalState(1,k).(fields{i}) + nnback*neighbor_terminalState_t(j,k).(fields{i});
            end
        end

    end

    for i = [10:13,18:21]
        neighbor_terminalState(1,k).(fields{i}) = neighbor_terminalState(1,k).(fields{i})/nfront;
    end
    for i = [14:17,22:25]
        neighbor_terminalState(1,k).(fields{i}) = neighbor_terminalState(1,k).(fields{i})/nback;
    end
end