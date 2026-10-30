"""Create T-Rex up to guide 4a, the thermal path through the extraction guides, and print it."""

from mcdgs.instruments.trex.trex import TRex


def main() -> None:
    """Create the instrument, without what follows guide 4a, and print its summary."""
    trex = TRex()
    trex.bender = None
    trex.guide_4b = None
    trex.sample_pos = None
    print(trex)


if __name__ == "__main__":
    main()
