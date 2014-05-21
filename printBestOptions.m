function bestOpt=printBestOptions(opt)
fprintf('%.10g ',[opt.wtEdet opt.wtEdyn opt.wtEexc opt.wtEper opt.wtEreg opt.wtEapp opt.wtEori opt.lambda]); 
bestOpt=[opt.wtEdet opt.wtEdyn opt.wtEexc opt.wtEper opt.wtEreg opt.wtEapp opt.wtEori opt.lambda];
fprintf('\n');