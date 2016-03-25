

  % save as kml. pol='h' or 'v' 

  function save_kml (kml_file_name, Tb, lat_filtered, lon_filtered, ...
                                         alt_filtered, roll_filtered, pol) 

    minBin=170;
    maxBin=300;
    % elevation angle is a function of the plane's roll ange
    elev = 40 + roll_filtered;
    % AGL = Altitude above ground level
    AGL = alt_filtered;
    % determine footprint size based on altitude
    slant = AGL./ cosd(elev);
    % beam width for radiometer specified in SLAP documentation
    if strcmp(pol, 'h') 
        beam_width = 18.8;
        poltype = '_h';
    else
        beam_width = 18.1;
        poltype = '_v';
    end
    %YDT footprintm = (slant*tand(beam_width/2))/cosd(elev);
    footprintm = (slant*tand(beam_width/2))./cosd(elev);

    r_earth = 6378e3;
    % convert footprint radius in meters to degrees longitude
    footprint_major = footprintm / (pi*r_earth) * 180;
    % the 6-deg averaged footprint, where each individual -3 dB footprint
    % is an ellipse, is approximated as a circle
    footprint_minor = footprint_major;
    %% generate KML file to plot on Google Earth with lat&lon values and Tb/SM data
    
    'here in Google Earth'
    %keyboard
    % create KML file with the name specified
    fida=fopen(kml_file_name,'w');
    
    line1=['<?xml version="1.0" encoding="UTF-8"?>'];
    line2=['<kml xmlns="http://earth.google.com/kml/2.0">'];
    line3=['<Document>'];
    line4=[strcat('<name>',kml_file_name,'</name>')];
    line6=['<Style id="air"><icon><href>root://icons/palette-4.png?x=160&amp;y=0&amp;w=32&amp;h=32</href></icon>'];
    line7=['<LineStyle><color>FF0000FF</color><width>2.0</width></LineStyle></Style>'];
    line8=['<Style id="gnd"><icon><href>root://icons/palette-4.png?x=160&amp;y=0&amp;w=32&amp;h=32</href></icon>'];
    line9=['<LineStyle><color>FFBBFF00</color><width>2.0</width></LineStyle></Style>'];
    line10=['<IconStyle><color>FFFFFFFF</color><scale>0.7</scale></IconStyle>'];
    line11=['<LabelStyle><color>FFFFFFFF</color><scale>0.7</scale></LabelStyle></Style>'];
    line_folder=['</Folder>'];
    line_last=['</Document></kml>'];
    
    
    cls=['FF0000AA';'FF0000D4';'FF0000FF';'FF002AFF';'FF0055FF';'FF0080FF';'FF00AAFF'; ...
        'FF00D4FF';'FF00FFFF';'FF2AFFD4';'FF54FFAA';'FF80FF80';'FFAAFF55';'FFD4FF2A'; ...
        'FFFFFF00';'FFFFD400';'FFFFAA00';'FFFF8000';'FFFF5500';'FFFF2A00';'FFFF0000';'50000000';'50000000'];

    cls =     flipud(cls);
    
    % to have more partitions of smaller ranges of values plotted on Google Earth
    %   cls =  ['50F00014'; '50F01414'; '50F02814'; '50F03C14'; '50F05014'; ...
    %      '50F06414'; '50F07814'; '50F08C14'; '50F0A014'; '50F0B414'; '50F0C814'; ...
    %      '50F0DC14'; '50F0F014'; '50F0FF14'; '50B4FF14'; '5078FF14'; '5014F000'; ...
    %      '5014F014'; '5014F028'; '5014F03C'; '5014F050'; '5014F064'; '5014F078'; ...
    %      '5014F08C'; '5014F0A0'; '5014F0B4'; '5014F0C8'; '5014F0DC'; '5014F0F0'; ...
    %      '5014F0FF'; '5014B4FF'; '501478FF'; '50143CFF'; '501400FF'; '501400FF'; ...
    %      '501400DC'; '501400C8'; '501400B4'; '501400A0'; '5014008C'; '50140078';'50140064'];
    
    fprintf(fida,'%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n', ...
        line1,line2,line3,line4,line6,line7,line8,line9);
    
    %-----------------------------------
    % bin into value ranges (20 colors and 2 outside of range)
    %-----------------------------------
    
    deltaB=(maxBin-minBin)./20;
    f0=find(Tb<minBin); % | isnan(Tb));
    f1=find(Tb>=minBin & Tb<minBin+deltaB*1);
    f2=find(Tb>=minBin+deltaB*1 & Tb<minBin+deltaB*2);
    f3=find(Tb>=minBin+deltaB*2 & Tb<minBin+deltaB*3);
    f4=find(Tb>=minBin+deltaB*3 & Tb<minBin+deltaB*4);
    f5=find(Tb>=minBin+deltaB*4 & Tb<minBin+deltaB*5);
    f6=find(Tb>=minBin+deltaB*5 & Tb<minBin+deltaB*6);
    f7=find(Tb>=minBin+deltaB*6 & Tb<minBin+deltaB*7);
    f8=find(Tb>=minBin+deltaB*7 & Tb<minBin+deltaB*8);
    f9=find(Tb>=minBin+deltaB*8 & Tb<=minBin+deltaB*9);
    f10=find(Tb>=minBin+deltaB*9 & Tb<=minBin+deltaB*10);
    f11=find(Tb>=minBin+deltaB*10 & Tb<=minBin+deltaB*11);
    f12=find(Tb>=minBin+deltaB*11 & Tb<=minBin+deltaB*12);
    f13=find(Tb>=minBin+deltaB*12 & Tb<=minBin+deltaB*13);
    f14=find(Tb>=minBin+deltaB*13 & Tb<=minBin+deltaB*14);
    f15=find(Tb>=minBin+deltaB*14 & Tb<=minBin+deltaB*15);
    f16=find(Tb>=minBin+deltaB*15 & Tb<=minBin+deltaB*16);
    f17=find(Tb>=minBin+deltaB*16 & Tb<=minBin+deltaB*17);
    f18=find(Tb>=minBin+deltaB*17 & Tb<=minBin+deltaB*18);
    f19=find(Tb>=minBin+deltaB*18 & Tb<=minBin+deltaB*19);
    f20=find(Tb>=minBin+deltaB*19 & Tb<=minBin+deltaB*20);
    f21=find(Tb>maxBin);
    % to have more partitions of smaller ranges of values plotted on Google Earth
    % f24=find(Tb>=minBin+deltaB*23 & Tb<=minBin+deltaB*24);
    % f25=find(Tb>=minBin+deltaB*24 & Tb<=minBin+deltaB*25);
    % f26=find(Tb>=minBin+deltaB*25 & Tb<=minBin+deltaB*26);
    % f27=find(Tb>=minBin+deltaB*26 & Tb<=minBin+deltaB*27);
    % f28=find(Tb>=minBin+deltaB*27 & Tb<=minBin+deltaB*28);
    % f29=find(Tb>=minBin+deltaB*28 & Tb<=minBin+deltaB*29);
    % f30=find(Tb>=minBin+deltaB*29 & Tb<=minBin+deltaB*30);
    % f31=find(Tb>=minBin+deltaB*30 & Tb<=minBin+deltaB*31);
    % f32=find(Tb>=minBin+deltaB*31 & Tb<=minBin+deltaB*32);
    % f33=find(Tb>=minBin+deltaB*32 & Tb<=minBin+deltaB*33);
    % f34=find(Tb>=minBin+deltaB*33 & Tb<=minBin+deltaB*34);
    % f35=find(Tb>=minBin+deltaB*34 & Tb<=minBin+deltaB*35);
    % f36=find(Tb>=minBin+deltaB*35 & Tb<=minBin+deltaB*36);
    % f37=find(Tb>=minBin+deltaB*36 & Tb<=minBin+deltaB*37);
    % f38=find(Tb>=minBin+deltaB*37 & Tb<=minBin+deltaB*38);
    % f39=find(Tb>=minBin+deltaB*38 & Tb<=minBin+deltaB*39);
    % f40=find(Tb>=minBin+deltaB*39 & Tb<=minBin+deltaB*40);
    % f41=find(Tb>=minBin+deltaB*40 & Tb<=minBin+deltaB*41);
    % f42=find(Tb>=minBin+deltaB*41 & Tb<=minBin+deltaB*42);
    % f43=find(Tb>maxBin);
    %
    
    %keyboard
    
    for j=1:22
        fprintf(fida,'%s %s %s %s %s\r','   <Placemark><name>TB Range: ',num2str((minBin+deltaB*(j-2))) ...
            ,'-',num2str((minBin+deltaB*(j-1))),'(K)</name>');
        fprintf(fida,'%s\r\n','       <Style>');
        fprintf(fida,'%s\r\n','           <LineStyle>');
        fprintf(fida,'%s %s %s\r\n','               <color>',cls(j,:),'</color>');
        fprintf(fida,'%s\r\n','               <width> 0 </width>');
        fprintf(fida,'%s\r\n','           </LineStyle>');
        fprintf(fida,'%s\r\n','           <PolyStyle>');
        fprintf(fida,'%s %s %s\r\n','               <color>',cls(j,:),'</color>');
        fprintf(fida,'%s\r\n','           </PolyStyle>');
        fprintf(fida,'%s\r\n','       </Style>');
        fprintf(fida,'%s\r\n','     <MultiGeometry>');
        fbin = eval(genvarname(['f',num2str(j-1)]));
        
        if ~isempty(fbin)
            for n=1:length(fbin)
                
                fprintf(fida,'%s\r\n','       <Polygon>');
                fprintf(fida,'%s\r\n','           <altitudeMode>relativetoground</altitudeMode>');
                fprintf(fida,'%s\r\n','           <outerBoundaryIs>');
                fprintf(fida,'%s\r\n','               <LinearRing>');
                fprintf(fida,'%s\r\n','                   <coordinates>');

                % create ellipse of observations
                for theta = 0:15:360
                    fprintf(fida,'%s%6f,%6f,%0f\r','                   ',lon_filtered(fbin(n))+footprint_minor(fbin(n))*cosd(theta),lat_filtered(fbin(n))+footprint_major(fbin(n))*sind(theta), 6000);
                end
                
                fprintf(fida,'%s\r\n','                   </coordinates>');
                fprintf(fida,'%s\r\n','               </LinearRing>');
                fprintf(fida,'%s\r\n','           </outerBoundaryIs>');
                fprintf(fida,'%s\r\n','       </Polygon>');
            end
        end
        
        fprintf(fida,'%s\r\n','     </MultiGeometry>');
        fprintf(fida,'%s\r\n','   </Placemark>');
        
    end
    
    fprintf(fida,'%s\r\n',line_last);
    fclose(fida);


end  % function 


