function ty_data = initialise_ty_data(R0, Fz0)


% Tyre structure data initialization

ty_data.FZ0             = Fz0; % Fz0  % Normal load
ty_data.R0              = R0; % R0  % nominal radius
ty_data.pCx1            = 1.5; %1; % pCx1
ty_data.pDx1            = 2.35; %1; % pDx1
ty_data.pDx2            = 0; % pDx2
ty_data.pDx3            = 0; % pDx3
ty_data.pEx1            = 0; % pEx1
ty_data.pEx2            = 0; % pEx2
ty_data.pEx3            = 0; % pEx3
ty_data.pEx4            = 0; % pEx4
ty_data.pHx1            = 0; % pHx1
ty_data.pHx2            = 0; % pHx2
ty_data.pKx1            = 50; %1; % pKx1
ty_data.pKx2            = 0; % pKx2
ty_data.pKx3            = 0; % pKx3
ty_data.pVx1            = 0; % pVx1
ty_data.pVx2            = 0; % pVx2
ty_data.Fz01            = 0; % Fz01
ty_data.pCy1            = 0; % pCy1
ty_data.pDy1            = 0; % pDy1
ty_data.pDy2            = 0; % pDy2
ty_data.pDy3            = 0; % pDy3
ty_data.pEy1            = 0; % pEy1
ty_data.pEy2            = 0; % pEy2
ty_data.pEy3            = 0; % pEy3
ty_data.pEy4            = 0; % pEy4
ty_data.pHy1            = 0; % pHy1
ty_data.pHy2            = 0; % pHy2
ty_data.pHy3            = 0; % pHy3
ty_data.pKy1            = 0; % pKy1
ty_data.pKy2            = 0; % pKy2
ty_data.pKy3            = 0; % pKy3
ty_data.pVy1            = 0; % pVy1
ty_data.pVy2            = 0; % pVy2
ty_data.pVy3            = 0; % pVy3
ty_data.pVy4            = 0; % pVy4
ty_data.qBz1            = 0; % qBz1
ty_data.qBz10           = 0; % qBz10
ty_data.qBz2            = 0; % qBz2
ty_data.qBz3            = 0; % qBz3
ty_data.qBz4            = 0; % qBz4
ty_data.qBz5            = 0; % qBz5
ty_data.qBz9            = 0; % qBz9
ty_data.qCz1            = 0; % qCz1
ty_data.qDz1            = 0; % qDz1
ty_data.qDz2            = 0; % qDz2
ty_data.qDz3            = 0; % qDz3
ty_data.qDz4            = 0; % qDz4
ty_data.qDz6            = 0; % qDz6
ty_data.qDz7            = 0; % qDz7
ty_data.qDz8            = 0; % qDz8
ty_data.qDz9            = 0; % qDz9
ty_data.qEz1            = 0; % qEz1
ty_data.qEz2            = 0; % qEz2
ty_data.qEz3            = 0; % qEz3
ty_data.qEz4            = 0; % qEz4
ty_data.qEz5            = 0; % qEz5
ty_data.qHz1            = 0; % qHz1
ty_data.qHz2            = 0; % qHz2
ty_data.qHz3            = 0; % qHz3
ty_data.qHz4            = 0; % qHz4
ty_data.rBx2            = 0; % rBx2
ty_data.rBy1            = 0; % rBy1
ty_data.rBy2            = 0; % rBy2
ty_data.rBy3            = 0; % rBy3
ty_data.rCx1            = 0; % rCx1
ty_data.rCy1            = 0; % rCy1
ty_data.rHx1            = 0; % rHx1
ty_data.rHy1            = 0; % rHy1
ty_data.rVy1            = 0; % rVy1
ty_data.rVy2            = 0; % rVy2
ty_data.rVy3            = 0; % rVy3
ty_data.rVy4            = 0; % rVy4
ty_data.rVy5            = 0; % rVy5
ty_data.rVy6            = 0; % rVy6



% scaling factor
ty_data.LCX             = 1; % LCX
ty_data.LCY             = 1; % LCY
ty_data.LEX             = 1; % LEX
ty_data.LEY             = 1; % LEY
ty_data.LFZ0            = 1; % LFZ0
ty_data.LGAMMAY         = 1; % LGAMMAY
ty_data.LHX             = 1; % LHX
ty_data.LHY             = 1; % LHY
ty_data.LKA             = 1; % LKA
ty_data.LKXK            = 1; % LKXK
ty_data.LMUX            = 1; % LMUX
ty_data.LMUY            = 1; % LMUY
ty_data.LVX             = 1; % LVX
ty_data.LVY             = 1; % LVY
ty_data.LVYK            = 1; % LVYK
ty_data.LXA             = 1; % LXA
ty_data.LYK             = 1; % LYK

end