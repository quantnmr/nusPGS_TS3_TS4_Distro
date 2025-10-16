# NUS Poisson Gap Sampling for TopSpin 3 & 4 v4.0

This software provides **Non-Uniform Sampling (NUS) Poisson Gap Sampling** for Bruker TopSpin 3.0+ and 4.0+ NMR spectrometers, reducing experimental time while maintaining spectral quality.

## Recommended Macros: PGS3 and PGS4

**We recommend using the simplified `PGS3` (TopSpin 3) or `PGS4` (TopSpin 4) macros** for most users. These provide a streamlined interface with sensible defaults while maintaining full functionality. The legacy `nusPGS_TS3` and `nusPGS_TS4` macros remain available for users requiring finer control of sampling parameters and are considered "battle tested" for production use.

**Bug reports and feedback**: Please send to scott.robson@northwestern.edu or scott@rezolytics.com

## Installation

### Step 1: Get the Software

If you have the repository URL, clone it:
```bash
git clone https://github.com/quantnmr/nusPGS_TS3_TS4_Distro.git
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

# Install recommended macro (choose based on your TopSpin version)
cp PGS3 /your/topspin/location/exp/stan/nmr/au/src/  # TopSpin 3 (RECOMMENDED)
cp PGS4 /your/topspin/location/exp/stan/nmr/au/src/  # TopSpin 4 (RECOMMENDED)

# Optional: Install legacy macros for advanced control
cp nusPGS_TS3 /your/topspin/location/exp/stan/nmr/au/src/  # TopSpin 3 (legacy)
cp nusPGS_TS4 /your/topspin/location/exp/stan/nmr/au/src/  # TopSpin 4 (legacy)
```

**For Windows:**
```cmd
# Install pre-compiled executable
copy poissonv3.eex C:\your\topspin\location\prog\bin\poissonv3.exe

# Install recommended macro (choose based on your TopSpin version)
copy PGS3 C:\your\topspin\location\exp\stan\nmr\au\src\  # TopSpin 3 (RECOMMENDED)
copy PGS4 C:\your\topspin\location\exp\stan\nmr\au\src\  # TopSpin 4 (RECOMMENDED)

# Optional: Install legacy macros for advanced control
copy nusPGS_TS3 C:\your\topspin\location\exp\stan\nmr\au\src\  # TopSpin 3 (legacy)
copy nusPGS_TS4 C:\your\topspin\location\exp\stan\nmr\au\src\  # TopSpin 4 (legacy)
```

## Usage

### Step 4: Test Installation

```bash
cd /your/topspin/location/prog/bin/
./poissonv3 2 0 2 818 0.001 64 128 0 0
```
Should output two columns of numbers.

### Step 5: Use in TopSpin

**Using PGS3 or PGS4 (Recommended - Simplified Interface):**

1. **Setup pulse program** normally in TopSpin
2. **Set acquisition parameter**: `FnTYPE` = `non-uniform_sampling`
3. **Run macro**: Type `PGS3` (TopSpin 3) or `PGS4` (TopSpin 4)
4. **Answer prompt**:
   - Number of Points: Enter desired number (see Quick Reference below for suggestions)
5. **Start acquisition**: `zg`

The macro uses optimized defaults: random seed (auto), sine portion (2), tolerance (0.00001), and shuffle (yes).

**Using nusPGS_TS3 or nusPGS_TS4 (Legacy - Advanced Control):**

1. **Setup pulse program** normally in TopSpin
2. **Set acquisition parameter**: `FnTYPE` = `non-uniform_sampling`
3. **Run macro**: Type `nusPGS_TS3` (TS3) or `nusPGS_TS4` (TS4)
4. **Answer prompts**:
   - Random Seed: `0` (or custom value)
   - Sine Portion: `2` (recommended)
   - Number of Points: Enter desired number
   - Tolerance: `0.01` or `0.001`
   - Shuffle: `1` (yes) or `0` (no)
5. **Start acquisition**: `zg`

## Quick Reference

- **2D experiments**: Use 25-50% sampling
- **3D experiments**: Use ~10% sampling
- **4D experiments**: Use 1-2% sampling

## Troubleshooting

- **Compilation fails**: Install `gcc` and math libraries
- **Permission denied**: Run `chmod +x poissonv3`
- **File not found**: Check TopSpin installation paths
- **Macro not found**: Ensure you copied the macro to the correct `au/src/` directory

## Feedback and Bug Reports

We encourage all users to transition to **PGS3** and **PGS4** for the best experience. Please report any issues, bugs, or suggestions to:
- **scott.robson@northwestern.edu**
- **scott@rezolytics.com**

Your feedback helps improve the software for the entire NMR community!
