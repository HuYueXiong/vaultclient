project "SPIRV-Cross"
    location "SPIRV-Cross"
	
	--Settings
	kind "StaticLib"

	language "C++"
	staticruntime "On"
	flags { "MultiProcessorCompile" }
	warnings "Off"

	--Files to include
	files { "SPIRV-Cross/*.h", "SPIRV-Cross/*.hpp", "SPIRV-Cross/*.cpp", "SPIRV-Cross/include/**", "SPIRV-Cross/external/spirv-tools/source/**", "SPIRV-Cross/external/spirv-tools/include/**" }
	files { "SPIRV-Cross.lua" }

	--This project includes
	includedirs { "SPIRV-Cross", "SPIRV-Cross/include", "SPIRV-Cross/external/spirv-tools/include", "SPIRV-Cross/external/spirv-tools/source", "SPIRV-Cross/external/spirv-tools", "SPIRV-Cross/external/spirv-tools/external/spirv-headers/include" }

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
