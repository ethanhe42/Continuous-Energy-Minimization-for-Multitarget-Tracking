function printProgressHeader(scenario,sceneInfo,opt)

if sceneInfo.gtAvailable
    if opt.track3d
        printMessage(3,'\n - INFO (%3d)- | -------------  E N E R G Y    V A L U E S  -------------- |||  ---------------- M E T R I C S (3D)---------------- |||  ----------------- M E T R I C S (2D)--------------- ',scenario);
    else
        printMessage(3,'\n - INFO (%3d)- | -------------  E N E R G Y    V A L U E S  -------------- |||  ----------------- M E T R I C S (2D)---------------- |||',scenario);
    end
    printMessage(3,'\n               |                                                           |||                                                      ||| ');
    
    if opt.track3d
        printMessage(3,'\n git| it|  time| Energy |  Edet  | Edyn | Eexc | Eapp | Eper | Ereg | Eori ||| MOTA  MOTP| GT  MT  ML|  FP   FN IDs  FM|  Rcll  Prcn||| MOTA  MOTP| GT  MT  ML|  FP   FN IDs  FM|  Rcll  Prcn\n');
    else
        printMessage(3,'\n git| it|  time| Energy |  Edet  | Edyn | Eexc | Eapp | Eper | Ereg | Eori ||| MOTA  MOTP| GT  MT  ML|  FP   FN IDs  FM|  Rcll  Prcn\n');
    end
else
    printMessage(3,'\n - INFO (%3d)- | ---------------  E N E R G Y    V A L U E S  ------------------ ',scenario);
    printMessage(3,'\n               |                                                           ');
    printMessage(3,'\n git| it|  time| Energy |  Edet  | Edyn | Eexc | Eapp | Eper | Ereg | Eori |\n');
end
end