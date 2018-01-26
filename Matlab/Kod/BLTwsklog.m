function Lcmax = BLTwsklog(G1,kp1,ki1,kp2,ki2,plot)
    %podanie modelu obiektu 2x2:
    s=tf('s');
    G11=G1(1,3);
    G12=G1(1,1);
    G21=G1(2,1);
    G22=G1(2,3);
    %regulatory PI:
    R1=kp1*(1+(ki1/s));
    R2=kp2*(1+(ki2/s));
    
    G11R1=series(G11,R1); %pierwszy sk³adnik sumy
    G22R2=series(G22,R2); %drugi sk³adnik sumy
    G11R1ss=ss(G11R1);
    G22R2ss=ss(G22R2);
    W1=parallel(G11R1ss,G22R2ss); %pierwsze sumowanie
    W1ss=ss(W1);
    GGRR1=series(G11R1ss,G22R2ss); %trzeci sk³adnik sumy
    W2=parallel(GGRR1,W1ss); %drugie sumowanie
    G12G21=series(G12,G21);
    R1R2=series(R1,R2);
    GGRR2=series(G12G21,R1R2);
    GGRR2m=series(GGRR2,-1); %czwarty sk³adnik sumy (suma z plusem)
    GGRR2mss=ss(GGRR2m);
    W=parallel(W2,GGRR2mss); %trzecie sumowanie (ostatnie)
    % ch-ka czêstotliwoœciowa:
    wektw=logspace(-4,-1,1000); % 100 punktów od 10^-4 do 10^-1
    Wfr=freqresp(W,wektw); %wektor ch-ki czêstotliwoœciowej W
    Wabs=abs(Wfr);
    W1=parallel(W,1); %1+W
    Wfr1=freqresp(W1,wektw); %wektor ch-ki czêstotliwoœciowej 1+W
    W1abs=abs(Wfr1);
    Lclogabs=20*(log10(Wabs)-log10(W1abs)); %wektor wartoœci funkcji wskaŸnika
    Lcmax=max(Lclogabs); %wartoœæ kryterium BLT
    if(plot == 1)
        semilogx(wektw, Lclogabs(1,:));
        grid;
    end
end