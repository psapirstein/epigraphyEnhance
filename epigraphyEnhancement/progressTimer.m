classdef progressTimer < handle
    properties (SetAccess = protected)
        delchars = 0;
        reportStep = 1;
    end
    properties (Hidden = true)
        N = 1;
        startTime;
    end
    methods
        function hd = progressTimer(Num,delchars)
            hd.startTime = clock;
            hd.delchars = delchars;
            hd.N = Num;
            hd.reportStep = floor(Num/250);
            if hd.reportStep < 1, hd.reportStep = 1; end
        end
        function update(hd,n)
            fprintf(repmat(char(8),1,hd.delchars));
            duration = etime(clock,hd.startTime);
            percProgress = n/hd.N;
            sRemaining = duration/percProgress-duration;
            if sRemaining < 100
                hd.delchars = fprintf('%d%% complete (%d sec remaining)',...
                    floor(100*percProgress),floor(sRemaining));
            elseif sRemaining < 300
                hd.delchars = fprintf('%d%% complete (%.1f min remaining)',...
                    floor(100*percProgress),floor(sRemaining/6)/10);
            else
                hd.delchars = fprintf('%d%% complete (%d min remaining)',...
                    floor(100*percProgress),floor(sRemaining/60));
            end
        end
        function done(hd)
            fprintf(repmat(char(8),1,hd.delchars));
        end
    end
end