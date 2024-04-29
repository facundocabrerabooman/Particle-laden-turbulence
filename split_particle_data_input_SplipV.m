clfc
% data to split
load('/Volumes/landau1/TrCer_analysis_paper#1/exports/ddt_pairs/particles/trajsf_TrCer_1000_02_ddt_particle.mat')

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong)==1);
tracklong = tracklong(Ine);

for i = 1:numel(tracklong)

    if ismembertol(tracklong(i).Tf(end), 0.847, 1e-2) 
        i
        tracklong1 = tracklong(1:i);
        tracklong2 = tracklong(i:end);
        break
    end

end

tracklong = tracklong1;
save(['/Volumes/landau1/TrCer_analysis_paper#1/exports/particles/ddt' filesep 'trajsf1_TrCer_1000_02_ddt_particle.mat'],'tracklong')

tracklong = tracklong2;
save(['/Volumes/landau1/TrCer_analysis_paper#1/exports/particles/ddt' filesep 'trajsf2_TrCer_1000_02_ddt_particle.mat'],'tracklong')