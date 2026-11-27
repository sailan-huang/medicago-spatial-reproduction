const fs = require("fs");
const path = require("path");
const {
  AlignmentType,
  BorderStyle,
  Document,
  ExternalHyperlink,
  Footer,
  Header,
  HeadingLevel,
  ImageRun,
  LevelFormat,
  PageBreak,
  PageNumber,
  Packer,
  Paragraph,
  ShadingType,
  Table,
  TableCell,
  TableRow,
  TextRun,
  WidthType,
} = require("docx");

const ROOT = path.resolve(__dirname, "..");
const REPORT_MD = path.join(ROOT, "docs", "analysis_report.md");
const METADATA_PATH = path.join(ROOT, "results", "tables", "run_metadata.json");
const ASSET_DIR = path.join(ROOT, "build", "docx_assets");
const OUTPUT = path.join(
  ROOT,
  "docs",
  "medicago_spatial_reproduction_report.docx",
);
const CONTENT_WIDTH = 9026;
const CJK_FONT = "Microsoft YaHei";
const FONT_SET = {
  ascii: "Arial",
  hAnsi: "Arial",
  eastAsia: CJK_FONT,
  cs: "Arial",
};
const tableBorder = { style: BorderStyle.SINGLE, size: 1, color: "C9D2DC" };
const tableBorders = {
  top: tableBorder,
  bottom: tableBorder,
  left: tableBorder,
  right: tableBorder,
};

const metadata = JSON.parse(fs.readFileSync(METADATA_PATH, "utf8"));
const reportLines = fs.readFileSync(REPORT_MD, "utf8").split(/\r?\n/);

function inlineRuns(text) {
  const cleaned = text
    .replace(/\[([^\]]+)\]\(([^)]+)\)/g, "$1")
    .replace(/\*\*/g, "")
    .replace(/\*([^*]+)\*/g, "$1");
  const parts = cleaned.split(/(`[^`]+`)/g).filter(Boolean);
  return parts.map(
    (part) =>
      new TextRun(
        part.startsWith("`") && part.endsWith("`")
          ? {
              text: part.slice(1, -1),
              font: "Consolas",
              size: 18,
              color: "C7254E",
            }
          : { text: part },
      ),
  );
}

function pngDimensions(file) {
  const buffer = fs.readFileSync(file);
  if (buffer.toString("ascii", 1, 4) !== "PNG") {
    throw new Error(`Not a PNG file: ${file}`);
  }
  return {
    width: buffer.readUInt32BE(16),
    height: buffer.readUInt32BE(20),
  };
}

function figureParagraph(stem, title) {
  const imagePath = path.join(ASSET_DIR, `${stem}.png`);
  if (!fs.existsSync(imagePath)) {
    throw new Error(`Missing DOCX preview image: ${imagePath}`);
  }
  const dimensions = pngDimensions(imagePath);
  const width = 600;
  const height = Math.round((width * dimensions.height) / dimensions.width);
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    keepNext: true,
    spacing: { before: 160, after: 80 },
    children: [
      new ImageRun({
        type: "png",
        data: fs.readFileSync(imagePath),
        transformation: { width, height },
        altText: {
          title,
          description: title,
          name: stem,
        },
      }),
    ],
  });
}

function tableFromRows(rows) {
  const columnCount = Math.max(...rows.map((row) => row.length));
  const columnWidth = Math.floor(CONTENT_WIDTH / columnCount);
  const widths = Array(columnCount).fill(columnWidth);
  widths[widths.length - 1] += CONTENT_WIDTH - columnWidth * columnCount;
  return new Table({
    width: { size: CONTENT_WIDTH, type: WidthType.DXA },
    columnWidths: widths,
    rows: rows.map(
      (row, rowIndex) =>
        new TableRow({
          tableHeader: rowIndex === 0,
          children: widths.map(
            (width, columnIndex) =>
              new TableCell({
                borders: tableBorders,
                width: { size: width, type: WidthType.DXA },
                margins: { top: 70, bottom: 70, left: 90, right: 90 },
                shading:
                  rowIndex === 0
                    ? { fill: "D9EAF7", type: ShadingType.CLEAR }
                    : undefined,
                verticalAlign: "center",
                children: [
                  new Paragraph({
                    spacing: { after: 0 },
                    children: [
                      new TextRun({
                        text: row[columnIndex] || "",
                        bold: rowIndex === 0,
                        size: 15,
                        font: FONT_SET,
                      }),
                    ],
                  }),
                ],
              }),
          ),
        }),
    ),
  });
}

function parseReport(lines) {
  const children = [];
  let paragraphBuffer = [];
  let captionBuffer = [];

  const flushParagraph = () => {
    if (!paragraphBuffer.length) return;
    children.push(
      new Paragraph({
        alignment: AlignmentType.JUSTIFIED,
        spacing: { after: 100, line: 276 },
        children: inlineRuns(paragraphBuffer.join(" ")),
      }),
    );
    paragraphBuffer = [];
  };

  const flushCaption = () => {
    if (!captionBuffer.length) return;
    children.push(
      new Paragraph({
        style: "Caption",
        alignment: AlignmentType.CENTER,
        spacing: { before: 40, after: 150 },
        children: inlineRuns(captionBuffer.join(" ")),
      }),
    );
    captionBuffer = [];
  };

  const flushAll = () => {
    flushParagraph();
    flushCaption();
  };

  for (let index = 0; index < lines.length; index += 1) {
    const raw = lines[index].trim();
    if (!raw) {
      flushAll();
      continue;
    }
    if (raw.startsWith("# ")) {
      continue;
    }
    if (raw.startsWith("## ")) {
      flushAll();
      children.push(
        new Paragraph({
          heading: HeadingLevel.HEADING_1,
          pageBreakBefore: ["结果", "与论文和作者代码的关系"].includes(
            raw.slice(3),
          ),
          children: [new TextRun(raw.slice(3))],
        }),
      );
      continue;
    }
    if (raw.startsWith("### ")) {
      flushAll();
      children.push(
        new Paragraph({
          heading: HeadingLevel.HEADING_2,
          children: [new TextRun(raw.slice(4))],
        }),
      );
      continue;
    }

    const figureMatch = raw.match(
      /^\[([^\]]+)\]\(\.\.\/results\/figures\/([^/()]+)\.pdf\)$/,
    );
    if (figureMatch) {
      flushAll();
      children.push(figureParagraph(figureMatch[2], figureMatch[1]));
      continue;
    }

    if (raw.startsWith("|")) {
      flushAll();
      const tableRows = [];
      while (index < lines.length && lines[index].trim().startsWith("|")) {
        const row = lines[index]
          .trim()
          .replace(/^\||\|$/g, "")
          .split(/(?<!\\)\|/)
          .map((cell) => cell.trim().replace(/\\\|/g, "|"));
        if (!row.every((cell) => /^:?-{3,}:?$/.test(cell))) {
          tableRows.push(row);
        }
        index += 1;
      }
      index -= 1;
      children.push(tableFromRows(tableRows));
      children.push(new Paragraph({ spacing: { after: 80 } }));
      continue;
    }

    const bulletMatch = raw.match(/^- (.+)$/);
    if (bulletMatch) {
      flushAll();
      children.push(
        new Paragraph({
          numbering: { reference: "bullets", level: 0 },
          spacing: { after: 60, line: 260 },
          children: inlineRuns(bulletMatch[1]),
        }),
      );
      continue;
    }

    const numberMatch = raw.match(/^\d+\. (.+)$/);
    if (numberMatch) {
      flushAll();
      children.push(
        new Paragraph({
          numbering: { reference: "numbers", level: 0 },
          spacing: { after: 60, line: 260 },
          children: inlineRuns(numberMatch[1]),
        }),
      );
      continue;
    }

    if (/^图\s+\d/.test(raw)) {
      flushParagraph();
      captionBuffer.push(raw);
      continue;
    }

    if (captionBuffer.length) {
      captionBuffer.push(raw);
      continue;
    }
    paragraphBuffer.push(raw);
  }
  flushAll();
  return children;
}

const titlePage = [
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 1400, after: 300 },
    children: [
      new TextRun({
        text: "苜蓿根瘤空间转录组",
        bold: true,
        size: 42,
        color: "17365D",
        font: FONT_SET,
      }),
    ],
  }),
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { after: 280 },
    children: [
      new TextRun({
        text: "单切片局部复现报告",
        bold: true,
        size: 34,
        color: "17365D",
        font: FONT_SET,
      }),
    ],
  }),
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { after: 100 },
    children: [
      new TextRun({
        text: "Spatiotemporal co-transcriptomics in Medicago nodules",
        italics: true,
        size: 22,
        color: "555555",
      }),
    ],
  }),
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 800, after: 80 },
    children: [new TextRun({ text: "作者：Sailan Huang", size: 22 })],
  }),
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { after: 80 },
    children: [new TextRun({ text: "日期：2026年11月27日", size: 22 })],
  }),
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { after: 80 },
    children: [
      new ExternalHyperlink({
        link: "https://github.com/sailan-huang/medicago-spatial-reproduction",
        children: [
          new TextRun({
            text: "github.com/sailan-huang/medicago-spatial-reproduction",
            color: "2F75B5",
            underline: {},
            size: 20,
          }),
        ],
      }),
    ],
  }),
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 500 },
    children: [
      new TextRun({
        text: `代码基线：${metadata.git_revision || "not available"}`,
        color: "777777",
        size: 18,
      }),
    ],
  }),
  new Paragraph({ children: [new PageBreak()] }),
];

const document = new Document({
  creator: "Sailan Huang",
  title: "苜蓿根瘤空间转录组单切片复现报告",
  subject: "Cell 2026 Medicago spatial transcriptomics reproduction",
  description:
    "Single-section reanalysis with quality control, clustering, marker mapping and host-rhizobium signals.",
  styles: {
    default: {
      document: {
        run: { font: FONT_SET, size: 21, color: "222222" },
        paragraph: { spacing: { line: 276 } },
      },
    },
    paragraphStyles: [
      {
        id: "Title",
        name: "Title",
        basedOn: "Normal",
        next: "Normal",
        quickFormat: true,
        run: { size: 42, bold: true, color: "17365D", font: FONT_SET },
      },
      {
        id: "Heading1",
        name: "Heading 1",
        basedOn: "Normal",
        next: "Normal",
        quickFormat: true,
        run: { size: 30, bold: true, color: "17365D", font: FONT_SET },
        paragraph: {
          spacing: { before: 260, after: 140 },
          outlineLevel: 0,
          keepNext: true,
        },
      },
      {
        id: "Heading2",
        name: "Heading 2",
        basedOn: "Normal",
        next: "Normal",
        quickFormat: true,
        run: { size: 25, bold: true, color: "2F5597", font: FONT_SET },
        paragraph: {
          spacing: { before: 200, after: 100 },
          outlineLevel: 1,
          keepNext: true,
        },
      },
      {
        id: "Caption",
        name: "Caption",
        basedOn: "Normal",
        next: "Normal",
        quickFormat: true,
        run: { size: 17, italics: true, color: "555555", font: FONT_SET },
        paragraph: { spacing: { after: 120 } },
      },
    ],
  },
  numbering: {
    config: [
      {
        reference: "bullets",
        levels: [
          {
            level: 0,
            format: LevelFormat.BULLET,
            text: "•",
            alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 560, hanging: 280 } } },
          },
        ],
      },
      {
        reference: "numbers",
        levels: [
          {
            level: 0,
            format: LevelFormat.DECIMAL,
            text: "%1.",
            alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 560, hanging: 280 } } },
          },
        ],
      },
    ],
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: 11906, height: 16838 },
          margin: { top: 1300, right: 1440, bottom: 1300, left: 1440 },
        },
      },
      headers: {
        default: new Header({
          children: [
            new Paragraph({
              alignment: AlignmentType.RIGHT,
              border: {
                bottom: {
                  style: BorderStyle.SINGLE,
                  size: 4,
                  color: "B4C6E7",
                  space: 4,
                },
              },
              children: [
                new TextRun({
                  text: "DCS｜苜蓿根瘤空间转录组复现",
                  size: 17,
                  color: "666666",
                  font: FONT_SET,
                }),
              ],
            }),
          ],
        }),
      },
      footers: {
        default: new Footer({
          children: [
            new Paragraph({
              alignment: AlignmentType.CENTER,
              children: [
                new TextRun({
                  text: "第 ",
                  size: 17,
                  color: "666666",
                  font: FONT_SET,
                }),
                new TextRun({
                  children: [PageNumber.CURRENT],
                  size: 17,
                  color: "666666",
                }),
                new TextRun({
                  text: " 页",
                  size: 17,
                  color: "666666",
                  font: FONT_SET,
                }),
              ],
            }),
          ],
        }),
      },
      children: [...titlePage, ...parseReport(reportLines)],
    },
  ],
});

Packer.toBuffer(document).then((buffer) => {
  fs.writeFileSync(OUTPUT, buffer);
  console.log(OUTPUT);
});

