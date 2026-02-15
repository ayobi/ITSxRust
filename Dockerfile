# -------- builder --------
FROM rust:1-bookworm AS builder

WORKDIR /app

COPY Cargo.toml Cargo.lock ./
COPY src ./src

RUN cargo build --release --locked

# -------- runtime --------
FROM debian:bookworm-slim AS runtime

RUN apt-get update \
  && apt-get install -y --no-install-recommends hmmer ca-certificates \
  && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/itsxrust /usr/local/bin/itsxrust

ENTRYPOINT ["itsxrust"]
CMD ["--help"]
