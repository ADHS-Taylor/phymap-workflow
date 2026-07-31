# setup.ps1
# Automates the setup of phymap-workflow and phymapr R package dependencies.

$ErrorActionPreference = "Stop"
Write-Host "=== Starting phymap-workflow Setup ===" -ForegroundColor Green

# 1. Locate Conda
$condaBat = ""
if (Test-Path "$env:USERPROFILE\miniconda3\condabin\conda.bat") {
    $condaBat = "$env:USERPROFILE\miniconda3\condabin\conda.bat"
} elseif (Test-Path "$env:USERPROFILE\anaconda3\condabin\conda.bat") {
    $condaBat = "$env:USERPROFILE\anaconda3\condabin\conda.bat"
} else {
    $condaCmd = Get-Command conda -ErrorAction SilentlyContinue
    if ($condaCmd) {
        $condaBat = $condaCmd.Source
    }
}

if (-not $condaBat) {
    Write-Error "Conda was not found in standard paths. Please ensure Miniconda/Anaconda is installed."
    Exit 1
}
Write-Host "[OK] Found Conda at: $condaBat" -ForegroundColor Green

# 2. Accept Conda Terms of Service (Default channels)
Write-Host "=== Accepting Conda Terms of Service ===" -ForegroundColor Cyan
& $condaBat tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
& $condaBat tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
& $condaBat tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2

# 3. Create/Update conda environment
Write-Host "=== Creating/Updating conda environment 'phymapr' ===" -ForegroundColor Cyan
& $condaBat create -c bioconda -c conda-forge -n phymapr snakemake-minimal -y
& $condaBat install -n phymapr -c conda-forge -c bioconda treetime pandas biopython -y

# 4. WSL Setup for RAxML-NG (Windows only)
if ($IsWindows -or $env:OS -like "*Windows*") {
    Write-Host "=== Checking WSL for RAxML-NG ===" -ForegroundColor Cyan
    $wslCheck = Get-Command wsl -ErrorAction SilentlyContinue
    if (-not $wslCheck) {
        Write-Warning "WSL was not found. If you are on Windows, RAxML-NG will not run unless installed manually."
    } else {
        # Check if raxml-ng is installed in WSL
        $raxmlCheck = wsl bash -c "which raxml-ng || test -f ~/.local/bin/raxml-ng && echo 'OK'"
        if ($raxmlCheck -like "*OK*") {
            Write-Host "[OK] RAxML-NG is already installed in WSL." -ForegroundColor Green
        } else {
            Write-Host "Installing RAxML-NG in WSL..." -ForegroundColor Yellow
            wsl bash -c "mkdir -p ~/.local/bin && curl -L -o /tmp/raxml-ng.zip https://github.com/amkozlov/raxml-ng/releases/download/2.0.2/raxml-ng_v2.0.2_linux_x86_64.zip && python3 -m zipfile -e /tmp/raxml-ng.zip /tmp/raxml-ng-bin && cp /tmp/raxml-ng-bin/raxml-ng ~/.local/bin/ && chmod +x ~/.local/bin/raxml-ng"
            Write-Host "[OK] RAxML-NG installed in WSL." -ForegroundColor Green
        }
    }
}

# 5. R Package Dependencies Setup
Write-Host "=== Setting up R Package Dependencies ===" -ForegroundColor Cyan
$rscriptPath = ""
if (Test-Path "C:\Program Files\R\R-4.2.3\bin\Rscript.exe") {
    $rscriptPath = "C:\Program Files\R\R-4.2.3\bin\Rscript.exe"
} else {
    $rCmd = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($rCmd) {
        $rscriptPath = $rCmd.Source
    }
}

if (-not $rscriptPath) {
    Write-Warning "Rscript was not found. Please install R (v4.2+) and add it to your PATH."
} else {
    Write-Host "[OK] Found Rscript at: $rscriptPath" -ForegroundColor Green
    
    # Install CRAN packages from PPM snapshot
    Write-Host "Installing CRAN dependencies..." -ForegroundColor Yellow
    & $rscriptPath -e "install.packages(c('remotes', 'BiocManager', 'ape', 'dplyr', 'gganimate', 'gifski', 'lubridate', 'maps', 'phytools', 'sf', 'tibble', 'tidygeocoder', 'tidyr', 'ggspatial'), repos='https://packagemanager.posit.co/cran/2023-10-31', type='binary')"
    
    # Install Bioconductor packages
    Write-Host "Installing Bioconductor dependencies..." -ForegroundColor Yellow
    & $rscriptPath -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2023-10-31')); BiocManager::install(c('ggtree', 'treeio'), update=FALSE, ask=FALSE)"
    
    # Install local phymapr package
    $phymaprPath = Resolve-Path "../phymapr" -ErrorAction SilentlyContinue
    if ($phymaprPath) {
        $phymaprPathStr = $phymaprPath.Path.Replace('\', '/')
        Write-Host "Installing local phymapr R package from: $phymaprPathStr" -ForegroundColor Yellow
        & $rscriptPath -e "remotes::install_local('$phymaprPathStr', upgrade='never', force=TRUE)"
        Write-Host "[OK] Installed phymapr R package." -ForegroundColor Green
    } else {
        Write-Error "Local phymapr package folder not found at ../phymapr."
    }
}

Write-Host "=== Setup Complete! ===" -ForegroundColor Green
Write-Host "To run the pipeline:" -ForegroundColor Yellow
Write-Host "1. Activate conda environment: conda activate phymapr" -ForegroundColor White
Write-Host "2. (Windows) Ensure R is in path: `$env:PATH = 'C:\Program Files\R\R-4.2.3\bin;' + `$env:PATH" -ForegroundColor White
Write-Host "3. Run Snakemake: snakemake -j 4" -ForegroundColor White
