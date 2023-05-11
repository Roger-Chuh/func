function testAB()
fromGlobal.alpha = -0.21;
fromGlobal.beta = 10;


toGlobal.alpha = 0.43;
toGlobal.beta = -50;

toRelative = calcRelative(fromGlobal, toGlobal);

glob = calcGlobal(fromGlobal, toRelative);

err = [glob.alpha - toGlobal.alpha glob.beta - toGlobal.beta]

end
function rel = calcRelative(fromGlobal, toGlobal)
	
		relAlpha = toGlobal.alpha - fromGlobal.alpha;
		relBeta = toGlobal.beta - fromGlobal.beta*exp(relAlpha);

        rel.alpha = relAlpha;
        rel.beta = relBeta;
		
	
end
function glob = calcGlobal(fromGlobal, toRelative)
	
	 toGlobalAlpha = toRelative.alpha + fromGlobal.alpha;
	toGlobalBeta = toRelative.beta + fromGlobal.beta*exp(toRelative.alpha);

        glob.alpha = toGlobalAlpha;
        glob.beta = toGlobalBeta;
		
	
end