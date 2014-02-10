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

	project "phylo"

		kind "ConsoleApp"
		language "C++"
		targetname "phylo"
		files { "src/**" }
		includedirs { "src" }
