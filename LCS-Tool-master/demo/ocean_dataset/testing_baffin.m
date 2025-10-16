%% testing baffin script read


% make figure and save

x = 0:100;
y = sin(x);
figure(1)
plot(x,y);
title_string = 'sine_wave';
title(title_string);

print(title_string,'-dpng');

% save data

save('sine_wave.mat','y');

exit