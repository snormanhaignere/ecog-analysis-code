StimOrder={};
for i=1:36
    if i<10
        StimOrder{i}=['ms1_1B_s0' num2str(i) '_cut.wav'];
    else
    StimOrder{i}=['ms1_1B_s' num2str(i) '_cut.wav'];
    end
end