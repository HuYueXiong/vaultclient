project "DirectXShaderCompiler"
    location "DirectXShaderCompiler"
	
	--Settings
	kind "StaticLib"

	language "C++"
	staticruntime "On"
	flags { "MultiProcessorCompile" }
	warnings "Off"

	--Files to include
	files { "DirectXShaderCompiler/lib/**", "DirectXShaderCompiler/include/**" }
	files { "DirectXShaderCompiler.lua" }

	--This project includes
	includedirs { "DirectXShaderCompiler/include", "DirectXShaderCompiler/lib/*" }

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
