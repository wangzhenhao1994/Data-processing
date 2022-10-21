intensity=containers.Map({1,2,3,4},{'1.9E+14','4.5E+14','8.9E+14','1.2E+15'});
data=containers.Map();
delay = containers.Map();
scanLength = 20*1000;
for i = [1,2,3]
    delay(intensity(i))=importdata(strcat('Delay',intensity(i),'.dat'));
end
for mass = [1,17,18]
    for i = [1,2,3]
        strcat('Time',"Mass",string(mass),'H2O','piezo',intensity(i),'.dat')
        data(strcat(string(mass),'_',intensity(i)))=importdata(strcat('Time',"Mass",string(mass),'H2O','piezo',intensity(i),'.dat'));
    end
end
d = delay(intensity(1));
1/(d(2)-d(1))
plot(d,data(strcat(string(18),'_',intensity(1))))

wcoherence(data(strcat(string(18),'_',intensity(1))),data(strcat(string(17),'_',intensity(1))),1/(d(2)-d(1))*1E15,'PhaseDisplayThreshold',0.7);