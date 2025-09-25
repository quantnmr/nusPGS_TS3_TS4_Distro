# NUS Poisson Gap Sampling for TopSpin 3 & 4

This software provides **Non-Uniform Sampling (NUS) Poisson Gap Sampling** for Bruker TopSpin 3.0+ and 4.0+ NMR spectrometers, reducing experimental time while maintaining spectral quality.

## Installation

### Step 1: Get the Software

If you have the repository URL, clone it:
```bash
git clone https://github.com/nomadiq/nusPGS_TS3_TS4_Distro.git
cd nusPGS_TS3_TS4_Distro
```

Or download and extract the files to a directory.

### Step 2: Compile the Executable (Linux/Unix only)

```bash
./cmppoiss
```

### Step 3: Install Files

**For Linux/Unix:**
```bash
# Install executable
cp poissonv3 /your/topspin/location/prog/bin/
chmod +x /your/topspin/location/prog/bin/poissonv3

# Install macro (choose one)
cp nusPGS_TS3 /your/topspin/location/exp/stan/nmr/au/src/  # TopSpin 3
cp nusPGS_TS4 /your/topspin/location/exp/stan/nmr/au/src/  # TopSpin 4
```

**For Windows:**
```cmd
# Install pre-compiled executable
copy poissonv3.eex C:\your\topspin\location\prog\bin\poissonv3.exe

# Install macro (choose one)
copy nusPGS_TS3 C:\your\topspin\location\exp\stan\nmr\au\src\  # TopSpin 3
copy nusPGS_TS4 C:\your\topspin\location\exp\stan\nmr\au\src\  # TopSpin 4
```

## Usage

### Step 4: Test Installation

```bash
cd /your/topspin/location/prog/bin/
./poissonv3 2 0 2 818 0.001 64 128 0 0
```
Should output two columns of numbers.

### Step 5: Use in TopSpin

1. **Setup pulse program** normally in TopSpin
2. **Set acquisition parameter**: `FnTYPE` = `non-uniform_sampling`
3. **Run macro**: Type `nusPGS_TS3` (TS3) or `nusPGS_TS4` (TS4)
4. **Answer prompts** with defaults:
   - Random Seed: `0`
   - Sine Portion: `2`
   - Number of Points: The default number inserted is around 10% of total (adjust as needed)
   - Tolerance: `0.001`
   - Shuffle: `1`
5. **Start acquisition**: `zg`

## Quick Reference

- **2D experiments**: Use 25-50% sampling
- **3D experiments**: Use ~10% sampling
- **4D experiments**: Use 1-2% sampling

## Troubleshooting

- **Compilation fails**: Install `gcc` and math libraries
- **Permission denied**: Run `chmod +x poissonv3`
- **File not found**: Check TopSpin installation paths

**Support**: scott@rezolytics.com
