solution "PhyloBlaster"

	configurations {"Release", "Debug"}
	
	location "build"
	targetdir "build"
	flags { "Symbols" }
	defines {"_CRT_SECURE_NO_WARNINGS"}

	configuration "Debug"
		defines {"DEBUG"}
	configuration "Release"
		flags {"Optimize"}
		
	project "phylolib"
		kind "StaticLib"
		language "C++"
		files { "src/phylo/**" }
		files { "src/mtw/**" }

		includedirs { "src/phylo" }
		includedirs { "src" }

	project "phylo"

		kind "ConsoleApp"
		language "C++"
		targetname "phylo"
		files { "src/main.cpp" }
		includedirs { "src/phylo", "src" }
		links { "phylolib" }

	project "treefactory"
		kind "ConsoleApp"
		language "C++"
		targetname "treefactory"
		includedirs { "src/phylo", "src" }
		files { "src/treefactory.cpp" }
		links { "phylolib" }

		