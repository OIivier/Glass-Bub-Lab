NEach = 2;
NFibr = 3;
NVert = 5;
NRadius = 3;

for iRadius = 1:NRadius
    for iFibr = 1:NFibr
        for iVert = 1:NVert
            for iEach=1:NEach
                Fibr = 0.2+0.2*(iFibr-1);
                Vert = (iVert-1)*4;
                MaxRadius = 1+1*(iRadius-1);
                MapName = MakeMap(Fibr,Vert,MaxRadius);
                SimFromMap(MapName);
            end
        end
    end
end
                
                