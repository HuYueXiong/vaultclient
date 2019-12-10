project "ShaderConductor"
    location "ShaderConductor"
	
	--Settings
	kind "StaticLib"

	language "C++"
	staticruntime "On"
	flags { "MultiProcessorCompile" }
	warnings "Off"

	--Files to include
	files { "ShaderConductor/Source/**", "ShaderConductor/Include/**" }
	files { "ShaderConductor.lua" }

	--This project includes
	includedirs { "ShaderConductor/Include", "ShaderConductor/Source/Core","ShaderConductor/Source/Tools","ShaderConductor/Source/Wrapper" }
    includedirs { "../3rdParty/DirectXShaderCompiler/include" }
	includedirs { "../3rdParty/SPIRV-Cross" }
	
	links { "SPIRV-Cross" }
	
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
