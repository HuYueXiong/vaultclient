project "HLSLcc"
    location "HLSLcc"
	
	--Settings
	kind "StaticLib"

	language "C++"
	staticruntime "On"
	flags { "MultiProcessorCompile" }
	warnings "Off"

	--Files to include
	files { "HLSLcc/src/**", "HLSLcc/include/**" }
	files { "HLSLcc.lua" }

	--This project includes
	includedirs { "HLSLcc", "HLSLcc/include", "HLSLcc/src/*" }

	symbols "On"
	
	-- filters
	filter { "configurations:Debug" }
		optimize "Debug"

	filter { "configurations:Debug", "system:Windows" }
		ignoredefaultlibraries { "libcmt" }

	filter { "configurations:Release" }
		optimize "Full"
		omitframepointer "On"
		flags { "NoBufferSecurityCheck" }

	targetdir "Lib"
	debugdir "Lib"
