# Atria

Atria is designed to trim adapters and low-quality bases of next-generation sequencing data. It infers the insert DNA precisely by integrating both adapter information and reverse-complementary properties of pair-end reads within a delicate decision tree.

If you use Atria, please cite the paper:
> Jiacheng Chuan, Aiguo Zhou, Lawrence Richard Hale, Miao He, Xiang Li, Atria: an ultra-fast and accurate trimmer for adapter and quality trimming, Gigabyte, 1, 2021  https://doi.org/10.46471/gigabyte.31

## Features

- Highly accurate Illumina adapter trimming
- Paired-end consensus calling
- Quality trimming
- Poly X tail trimming
- Hard clip 3' and 5' ends
- N tail trimming
- Filtering reads by the number of N bases
- Filtering reads by length
- Filtering reads by read complexity

- Multi-threading and multi-core parallel computing

## Contents

1. Installation guide

   1.1 [Release installation guide](docs/1.1.Release_installation_guide.md)

   1.2 [Install from source](docs/1.2.Install_from_source.md)

2. **[Atria trimming methods and usages](docs/2.Atria_trimming_methods_and_usages.md)**

3. [Benchmark toolkit](docs/3.Benchmark_toolkit.md)

4. [Atria development notes](docs/4.Development_notes.md)

5. **[Accuracy and speed benchmark](docs/5.Accuracy_and_speed_benchmark.md)**
