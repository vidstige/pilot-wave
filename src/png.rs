//! Minimal PNG writer: 8-bit RGB, uncompressed deflate blocks.
use std::io;
use std::path::Path;

use crate::image::Image;

const BLOCK: usize = 65535;

fn crc32(data: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for &byte in data {
        crc ^= byte as u32;
        for _ in 0..8 {
            crc = if crc & 1 == 1 { 0xEDB8_8320 ^ (crc >> 1) } else { crc >> 1 };
        }
    }
    !crc
}

fn adler32(data: &[u8]) -> u32 {
    let (mut a, mut b) = (1u32, 0u32);
    for &byte in data {
        a = (a + byte as u32) % 65521;
        b = (b + a) % 65521;
    }
    (b << 16) | a
}

fn chunk(out: &mut Vec<u8>, kind: &[u8; 4], data: &[u8]) {
    out.extend((data.len() as u32).to_be_bytes());
    let mut body = kind.to_vec();
    body.extend(data);
    out.extend(crc32(&body).to_be_bytes());
    out.splice(out.len() - 4..out.len() - 4, body);
}

/// Zlib stream made of stored (uncompressed) deflate blocks.
fn zlib_stored(raw: &[u8]) -> Vec<u8> {
    let mut out = vec![0x78, 0x01];
    let blocks = raw.chunks(BLOCK).collect::<Vec<_>>();
    for (i, block) in blocks.iter().enumerate() {
        let last = i + 1 == blocks.len();
        out.push(last as u8);
        let n = block.len() as u16;
        out.extend(n.to_le_bytes());
        out.extend((!n).to_le_bytes());
        out.extend(*block);
    }
    out.extend(adler32(raw).to_be_bytes());
    out
}

pub fn encode(image: &Image) -> Vec<u8> {
    let mut raw = Vec::with_capacity(image.height * (image.width * 3 + 1));
    for row in image.pixels.chunks(image.width) {
        raw.push(0);
        raw.extend(row.iter().flatten());
    }
    let mut header = Vec::new();
    header.extend((image.width as u32).to_be_bytes());
    header.extend((image.height as u32).to_be_bytes());
    header.extend([8, 2, 0, 0, 0]);
    let mut out = vec![0x89, b'P', b'N', b'G', 0x0D, 0x0A, 0x1A, 0x0A];
    chunk(&mut out, b"IHDR", &header);
    chunk(&mut out, b"IDAT", &zlib_stored(&raw));
    chunk(&mut out, b"IEND", &[]);
    out
}

pub fn write(path: impl AsRef<Path>, image: &Image) -> io::Result<()> {
    std::fs::write(path, encode(image))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc_matches_the_png_reference_value() {
        assert_eq!(crc32(b"IEND"), 0xAE42_6082);
    }

    #[test]
    fn adler_of_empty_stream_is_one() {
        assert_eq!(adler32(&[]), 1);
    }

    #[test]
    fn encoded_image_has_signature_and_end_chunk() {
        let bytes = encode(&Image::new(3, 2, [1, 2, 3]));
        assert_eq!(&bytes[..8], &[0x89, b'P', b'N', b'G', 0x0D, 0x0A, 0x1A, 0x0A]);
        assert_eq!(&bytes[bytes.len() - 8..bytes.len() - 4], b"IEND");
    }
}
