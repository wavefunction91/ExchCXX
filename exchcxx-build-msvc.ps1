# exchcxx-build-msvc.ps1
# Build and test ExchCXX on Windows with MSVC native cl.exe.
# Run from the ExchCXX repo root in a PowerShell.
#
# Usage:
#   .\exchcxx-build-msvc.ps1                  # Full build + test
#   .\exchcxx-build-msvc.ps1 -SkipConfigure   # Incremental build (skip cmake configure)
#   .\exchcxx-build-msvc.ps1 -SkipTests       # Build only, no tests

param(
    [switch]$SkipConfigure,
    [switch]$SkipTests
)

$ErrorActionPreference = "Stop"
$RepoRoot = $PSScriptRoot
$BuildDir = "$RepoRoot\build-msvc"
$QdkRoot = "C:\Users\v-lercole\src\qdk-chemistry"
$VcpkgInstalledDir = "$QdkRoot\vcpkg_installed"
$VcpkgTriplet = "x64-windows-static-md"

Write-Host "============================================" -ForegroundColor Cyan
Write-Host "  ExchCXX - Windows Build (MSVC)            " -ForegroundColor Cyan
Write-Host "============================================" -ForegroundColor Cyan
Write-Host "Repo root: $RepoRoot"
Write-Host ""

function Assert-Command($Name) {
    if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
        Write-Error "$Name not found in PATH. Please install it first."
        exit 1
    }
}

# --------------------------------------------------------------------------
# Set up MSVC environment
# --------------------------------------------------------------------------
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vsPath = & $vswhere -latest -products * -property installationPath 2>$null

# Set up MSVC developer environment (vcvarsall)
$vcvarsall = "$vsPath\VC\Auxiliary\Build\vcvarsall.bat"
if (Test-Path $vcvarsall) {
    $tempFile = [System.IO.Path]::GetTempFileName()
    cmd /c "`"$vcvarsall`" x64 && set > `"$tempFile`""
    Get-Content $tempFile | ForEach-Object {
        if ($_ -match "^([^=]+)=(.*)$") {
            [System.Environment]::SetEnvironmentVariable($matches[1], $matches[2], "Process")
        }
    }
    Remove-Item $tempFile
}

# Set up vcpkg environment (reuse qdk-chemistry's vcpkg_installed)
$vcpkgRoot = $null
$vsVcpkg = "$vsPath\VC\vcpkg\vcpkg.exe"
if (Test-Path $vsVcpkg) {
    $vcpkgRoot = Split-Path $vsVcpkg
} elseif ($env:VCPKG_INSTALLATION_ROOT -and (Test-Path "$($env:VCPKG_INSTALLATION_ROOT)\vcpkg.exe")) {
    $vcpkgRoot = $env:VCPKG_INSTALLATION_ROOT
} else {
    $vcpkgRoot = "$QdkRoot\vcpkg-tool"
}

$toolchainFile = "$vcpkgRoot\scripts\buildsystems\vcpkg.cmake"
$env:CMAKE_TOOLCHAIN_FILE = $toolchainFile
$env:VCPKG_TARGET_TRIPLET = $VcpkgTriplet
$env:VCPKG_INSTALLED_DIR = $VcpkgInstalledDir
$env:CMAKE_PREFIX_PATH = "$VcpkgInstalledDir\$VcpkgTriplet"
$env:PATH = "$VcpkgInstalledDir\$VcpkgTriplet\bin;$VcpkgInstalledDir\$VcpkgTriplet\debug\bin;$env:PATH"

# --------------------------------------------------------------------------
# Verify environment
# --------------------------------------------------------------------------
foreach ($cmd in @("cl", "cmake", "ninja")) {
    Assert-Command $cmd
}

if (-not $env:CMAKE_TOOLCHAIN_FILE) {
    Write-Error "CMAKE_TOOLCHAIN_FILE not set. Run qdk-chemistry's test-windows-build-msvc.ps1 -SkipCpp -SkipPython first to set up the environment."
    exit 1
}

Write-Host "=== Environment ===" -ForegroundColor Yellow
Write-Host "  cl:              $(Get-Command cl | Select-Object -ExpandProperty Source)"
Write-Host "  cmake:           $(Get-Command cmake | Select-Object -ExpandProperty Source)"
Write-Host "  ninja:           $(Get-Command ninja | Select-Object -ExpandProperty Source)"
Write-Host "  TOOLCHAIN_FILE:  $env:CMAKE_TOOLCHAIN_FILE"
Write-Host "  VCPKG_INSTALLED: $VcpkgInstalledDir"
Write-Host ""

# --------------------------------------------------------------------------
# Configure
# --------------------------------------------------------------------------

if (-not $SkipConfigure) {
    Write-Host "=== Configure ===" -ForegroundColor Yellow

    cmake -S "$RepoRoot" -B "$BuildDir" `
        -GNinja `
        -DCMAKE_POLICY_VERSION_MINIMUM="3.5" `
        -DFETCHCONTENT_QUIET=OFF `
        -DBUILD_SHARED_LIBS=OFF `
        -DCMAKE_BUILD_TYPE=Release `
        -DCMAKE_CXX_STANDARD=20 `
        -DCMAKE_CXX_STANDARD_REQUIRED=ON `
        -DCMAKE_C_COMPILER=cl `
        -DCMAKE_CXX_COMPILER=cl `
        -DCMAKE_CXX_FLAGS="/W3 /EHsc /Zc:preprocessor /wd4800 /wd4514 /wd4710 /wd4711 /wd4061 /wd4266 /wd4702 /w14018 /w14055 /w14100 /w14102 /w14127 /w14146 /w14242 /w14244 /w14245 /w14254 /w14267 /w14302 /w14306 /w14308 /w14310 /w14389 /w14509 /w14510 /w14512 /w14532 /w14533 /w14610 /w14611 /w14700 /w14701 /w14703 /w14789 /w14995 /w14996" `
        -DCMAKE_CXX_FLAGS_RELEASE="/O2 /DNDEBUG" `
        -DCMAKE_TOOLCHAIN_FILE="$env:CMAKE_TOOLCHAIN_FILE" `
        -DVCPKG_TARGET_TRIPLET=$VcpkgTriplet `
        -DVCPKG_INSTALLED_DIR="$VcpkgInstalledDir" `
        -DEXCHCXX_ENABLE_LIBXC=ON `
        -DEXCHCXX_ENABLE_CUDA=OFF `
        -DEXCHCXX_ENABLE_HIP=OFF `
        -DEXCHCXX_ENABLE_SYCL=OFF `
        -DEXCHCXX_ENABLE_TESTS=ON
    if ($LASTEXITCODE -ne 0) { Write-Error "CMake configure failed"; exit 1 }
    Write-Host "Configure succeeded." -ForegroundColor Green
} else {
    Write-Host "=== Skipping configure (incremental build) ===" -ForegroundColor DarkGray
    if (-not (Test-Path "$BuildDir\build.ninja")) {
        Write-Error "No existing build found at $BuildDir. Run without -SkipConfigure first."
        exit 1
    }
}

# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------
Write-Host ""
Write-Host "=== Build ===" -ForegroundColor Yellow
cmake --build "$BuildDir" --parallel 6 2>&1 *> "$BuildDir\build.log"
if ($LASTEXITCODE -ne 0) { Write-Error "Build failed (see $BuildDir\build.log)"; exit 1 }
Write-Host "Build succeeded." -ForegroundColor Green

# --------------------------------------------------------------------------
# Test
# --------------------------------------------------------------------------
if (-not $SkipTests) {
    Write-Host ""
    Write-Host "=== Test ===" -ForegroundColor Yellow
    $env:OMP_NUM_THREADS = 4
    Push-Location "$BuildDir"
    ctest --output-on-failure --verbose
    $ctestExit = $LASTEXITCODE
    Pop-Location
    if ($ctestExit -ne 0) {
        Write-Warning "Some tests failed (exit code: $ctestExit)"
    } else {
        Write-Host "All tests passed." -ForegroundColor Green
    }
}

Write-Host ""
Write-Host "============================================" -ForegroundColor Cyan
Write-Host "  Done!                                     " -ForegroundColor Cyan
Write-Host "============================================" -ForegroundColor Cyan
