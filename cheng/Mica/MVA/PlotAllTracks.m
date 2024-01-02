function PlotAllTracks(velm,decim)
L=numel(velm.good);

figure;hold on;CC=jet(L);arrayfun(@(X,I)(plot([1:numel(X.velf)]/32768/3.8e-3,X.velf,'Color',CC(I,:))),velm.data(velm.good(1:decim:L)),1:decim:L);

box on;
grid on;

% Create xlabel
xlabel('t/\tau_\eta','FontWeight','bold','FontSize',18);

% Create ylabel
ylabel('v [a.u.]','FontWeight','bold','FontSize',18);

% Create textbox
annotation(gcf,'textbox','String',{'\Gamma= 65; \Phi = 16.5'},...
	'FontSize',18,...
	'FitBoxToText','off',...
	'BackgroundColor',[0 0.749 0.749],...
	'Position',[0.6837 0.8376 0.21 0.06824]);

