VAULTSDK_HOME = path.translate(os.getenv('VAULTSDK_HOME') or '')

project "EuclideonVaultClient"
	-- Settings
	kind "packaging"
	system "android"

	files { "../../builds/**.json", "../../builds/**.md", "../../builds/**.otf", "../../builds/**.png", "../../builds/**.jpg" }
	files { "../../3rdParty/SDL2-2.0.8/lib/android/%{cfg.platform}/**.so" }
	files { "res/**.xml", "res/**.png" }
	files { "src/**.java"}

	vpaths {
		["libs/arm64-v8a/*"] = { "../../3rdParty/SDL2-2.0.8/lib/android/arm64/**" },
		["libs/x86_64/*"] = { "../../3rdParty/SDL2-2.0.8/lib/android/x64/**" },
		["assets/*"] = { "../../builds/**" },
		[""] = { "project.lua" },
	}

	files { "AndroidManifest.xml" }
	files { "build.xml" }
	files { "project.properties" }

	libdirs { "../../3rdParty/SDL2-2.0.8/lib/android/arm64" }
	links { "vaultClient", "SDL2", "main" }

	filter { "options:force-vaultsdk" }
		files { "../../../vault/buildscripts/android/src/**" }
		vpaths {
			["src/*"] = { "../../../vault/buildscripts/android/src/**" }
		}
	filter { "options:not force-vaultsdk" }
		files { "packaging.lua" }
		files { "%{VAULTSDK_HOME}/lib/android_arm64/src/**" }
		files { "%{VAULTSDK_HOME}/lib/android_arm64/libvaultSDK.so", "%{VAULTSDK_HOME}/lib/android_x64/libvaultSDK.so" }
		vpaths {
			["src/*"] = { "%{VAULTSDK_HOME}/lib/android_arm64/src/**", "src/**" },
			["libs/arm64-v8a/*"] = { "%{VAULTSDK_HOME}/lib/android_arm64/**.so" },
			["libs/x86_64/*"] = { "%{VAULTSDK_HOME}/lib/android_x64/**.so" },
		}
	filter {}

	dofile "../../3rdParty/udcore/bin/premake-bin/common-proj.lua"
