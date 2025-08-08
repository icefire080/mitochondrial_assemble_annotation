#!/usr/bin/env python3
import argparse
import os
import subprocess
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from collections import defaultdict
import uuid
import datetime
import shutil

class MitochondrialAnnotation:
    def __init__(self, query_file, reference_file, reference_annotation, output_dir, genetic_code=2, prefix=""):
        """
        初始化注释器
        Args:
            query_file: 待注释的目标序列FASTA文件
            reference_file: 参考序列FASTA文件
            reference_annotation: 参考序列GFF注释文件
            output_dir: 输出目录
            genetic_code: 遗传密码表编号
            prefix: 输出文件的前缀
        """
        self.query_file = query_file
        self.reference_file = reference_file
        self.reference_annotation = reference_annotation
        self.output_dir = output_dir
        self.genetic_code = genetic_code
        self.prefix = prefix
        
        # 读取序列
        self.query_seq = SeqIO.read(query_file, "fasta")
        self.reference_seq = SeqIO.read(reference_file, "fasta")
        
        # 若指定prefix,覆盖原序列ID
        if self.prefix:
            self.query_seq.id = self.prefix
        
        # 读取参考注释
        self.ref_features = self._parse_gff(reference_annotation)
        
        # 设置临时文件目录
        self.tmp_dir = os.path.join(output_dir, "tmp")
        os.makedirs(self.tmp_dir, exist_ok=True)
        
        # 设置注释时间戳和唯一标识符
        self.timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        self.uuid_base = str(uuid.uuid4())[:8]

    def _parse_gff(self, gff_file):
        """解析GFF文件并合并相同ID的CDS片段"""
        features = defaultdict(list)
        cds_groups = defaultdict(list)  # 按ID分组存储CDS
        
        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                    
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                phase = parts[7]
                
                attributes = {}
                for attr in parts[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value
                
                # 特殊处理CDS：按ID分组
                if feature_type == 'CDS':
                    feature_id = attributes.get('ID', attributes.get('Parent', None))
                    if not feature_id:
                        # 如果没有ID或Parent，作为独立特征处理
                        features[feature_type].append({
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'phase': phase,
                            'attributes': attributes,
                            'source': parts[1],
                            'is_segment': False
                        })
                    else:
                        # 存储分段信息
                        cds_groups[feature_id].append({
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'phase': phase,
                            'attributes': attributes,
                            'source': parts[1]
                        })
                else:
                    features[feature_type].append({
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'phase': phase,
                        'attributes': attributes,
                        'source': parts[1],
                        'is_segment': False
                    })
        
        # 合并分段CDS
        for feature_id, segments in cds_groups.items():
            if not segments:
                continue
                
            # 确保所有片段在同一链上
            strand = segments[0]['strand']
            for seg in segments[1:]:
                if seg['strand'] != strand:
                    print(f"Warning: CDS {feature_id} has segments on different strands. Using first strand.")
                    break
            
            # 按起始位置排序
            sorted_segments = sorted(segments, key=lambda x: x['start'])
            
            # 创建合并后的特征
            features['CDS'].append({
                'start': min(s['start'] for s in sorted_segments),
                'end': max(s['end'] for s in sorted_segments),
                'strand': strand,
                'phase': sorted_segments[0]['phase'],
                'attributes': sorted_segments[0]['attributes'].copy(),
                'source': sorted_segments[0]['source'],
                'segments': sorted_segments,  # 保存原始片段信息
                'is_segment': True
            })
        
        return features

    def run_whole_alignment(self):
        """
        执行全序列比对（使用mafft）
        返回比对后的参考序列和目标序列
        """
        # 准备输入文件（确保ID唯一）
        ref_input = os.path.join(self.tmp_dir, "ref_input.fasta")
        query_input = os.path.join(self.tmp_dir, "query_input.fasta")
        
        with open(ref_input, 'w') as f:
            f.write(f">ref\n{self.reference_seq.seq}\n")
        with open(query_input, 'w') as f:
            f.write(f">query\n{self.query_seq.seq}\n")
        # 创建合并的输入文件
        combined_input = os.path.join(self.tmp_dir, "combined_input.fasta")
    
        with open(combined_input, 'w') as f:
             f.write(f">ref\n{self.reference_seq.seq}\n")
             f.write(f">query\n{self.query_seq.seq}\n")
             
        # 运行mafft比对
        aligned_file = os.path.join(self.tmp_dir, "aligned.fasta")
        mafft_cmd = f"mafft --auto --reorder --anysymbol {combined_input}  > {aligned_file}"
        try:
            subprocess.run(mafft_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MAFFT alignment failed: {str(e)}")
        
        # 读取比对结果
        aligned_seqs = list(SeqIO.parse(aligned_file, "fasta"))
        ref_aligned = next(s for s in aligned_seqs if s.id == "ref")
        query_aligned = next(s for s in aligned_seqs if s.id == "query")
        
        return ref_aligned, query_aligned

    def _get_aligned_position(self, aligned_seq, orig_pos):
        """
        将原始序列位置转换为比对序列中的位置（考虑gap）
        """
        pos = 0
        aligned_pos = 0
        for base in aligned_seq.seq:
            if pos >= orig_pos:
                break
            if base != '-':
                pos += 1
            aligned_pos += 1
        return max(0, aligned_pos - 1)

    def extract_aligned_features(self, ref_aligned, query_aligned):
        """根据比对结果提取特征区域，处理分段CDS和负链"""
        features = defaultdict(list)
        
        for feature_type, feature_list in self.ref_features.items():
            for feature in feature_list:
                try:
                    # 处理分段CDS
                    if feature_type == 'CDS' and feature.get('is_segment', False):
                        # 分别提取每个片段
                        query_segments = []
                        ref_segments = []
                        segment_positions = []
                        
                        for segment in feature['segments']:
                            ref_start = self._get_aligned_position(ref_aligned, segment['start'])
                            ref_end = self._get_aligned_position(ref_aligned, segment['end'])
                            
                            # 提取序列并转换为字符串
                            query_seg_seq = str(query_aligned.seq[ref_start:ref_end+1])  # 修改点
                            ref_seg_seq = str(ref_aligned.seq[ref_start:ref_end+1])       # 修改点
                            
                            # 负链片段立即取反向互补
                            if feature['strand'] == '-':
                                query_seg_seq = str(Seq(query_seg_seq).reverse_complement())
                                ref_seg_seq = str(Seq(ref_seg_seq).reverse_complement())
                            
                            query_segments.append(query_seg_seq)
                            ref_segments.append(ref_seg_seq)
                            segment_positions.append((ref_start, ref_end))
                        
                        # 合并所有片段
                        query_seq = ''.join(query_segments)
                        ref_seq = ''.join(ref_segments)
                        
                        # 保存整个CDS特征
                        features[feature_type].append({
                            'start': min(p[0] for p in segment_positions) + 1,  # 1-based
                            'end': max(p[1] for p in segment_positions) + 1,
                            'strand': feature['strand'],
                            'phase': feature['phase'],
                            'attributes': feature['attributes'].copy(),
                            'query_seq': query_seq,
                            'ref_seq': ref_seq,
                            'segment_positions': segment_positions,
                            'segments': feature['segments']
                        })
                    else:
                        # 处理非分段特征
                        ref_start = self._get_aligned_position(ref_aligned, feature['start'])
                        ref_end = self._get_aligned_position(ref_aligned, feature['end'])
                        
                        # 提取序列并转换为字符串
                        query_feature_seq = str(query_aligned.seq[ref_start:ref_end+1])  # 修改点
                        ref_feature_seq = str(ref_aligned.seq[ref_start:ref_end+1])     # 修改点
                        
                        # 负链特征立即取反向互补
                        if feature['strand'] == '-':
                            query_feature_seq = str(Seq(query_feature_seq).reverse_complement())
                            ref_feature_seq = str(Seq(ref_feature_seq).reverse_complement())
                        
                        features[feature_type].append({
                            'start': ref_start + 1,
                            'end': ref_end + 1,
                            'strand': feature['strand'],
                            'phase': feature['phase'],
                            'attributes': feature['attributes'].copy(),
                            'query_seq': query_feature_seq,
                            'ref_seq': ref_feature_seq,
                            'segment_positions': [(ref_start, ref_end)]
                        })
                except Exception as e:
                    print(f"Warning: Failed to process feature {feature}: {str(e)}")
        
        return features

    def process_cds_features(self, features):
        """
        处理CDS特征：智能填充gap保证长度是3的倍数
        在参考序列的gap位置添加gap，使相连的gap为3的倍数
        """
        processed_features = []
        
        for feature in features.get('CDS', []):
            try:
                query_seq = feature['query_seq']
                ref_seq = feature['ref_seq']
                orig_len = len(query_seq)
                
                # 保存原始比对序列
                feature['original_query_seq'] = query_seq
                feature['original_ref_seq'] = ref_seq
                
                # 检查并填充gap
                if len(query_seq) % 3 != 0:
                    # 1. 识别参考序列中的gap区域
                    gap_regions = []
                    current_gap = None
                    for i, base in enumerate(ref_seq):
                        if base == '-':
                            if current_gap is None:
                                current_gap = [i, i]  # 开始新的gap区域
                            else:
                                current_gap[1] = i    # 扩展当前gap区域
                        else:
                            if current_gap is not None:
                                gap_regions.append(tuple(current_gap))
                                current_gap = None
                    
                    # 添加最后一个gap区域
                    if current_gap is not None:
                        gap_regions.append(tuple(current_gap))
                    
                    # 2. 计算需要填充的总gap数量
                    remainder = len(query_seq) % 3
                    padding_needed = 3 - remainder if remainder > 0 else 0
                    
                    # 3. 智能填充策略：
                    #    - 优先在已有gap区域末尾添加
                    #    - 如果无gap区域，在序列末尾添加
                    if gap_regions:
                        # 在最后一个gap区域后添加gap
                        last_gap_end = gap_regions[-1][1] + 1
                        padding_pos = last_gap_end
                    else:
                        # 无gap区域，在序列末尾添加
                        padding_pos = len(ref_seq)
                    
                    # 4. 在指定位置添加gap
                    padding = '-' * padding_needed
                    
                    # 在query序列插入位置添加N（避免破坏密码子）
                    new_query_seq = query_seq[:padding_pos] + 'N' * padding_needed + query_seq[padding_pos:]
                    
                    # 在ref序列相同位置添加gap
                    new_ref_seq = ref_seq[:padding_pos] + padding + ref_seq[padding_pos:]
                    
                    print(f"Adjusted CDS length: added {padding_needed} gaps at position {padding_pos} "
                        f"(original: {orig_len}, new: {len(new_query_seq)})")
                    
                    query_seq = new_query_seq
                    ref_seq = new_ref_seq
                
                # 保存处理后的序列
                feature['processed_query_seq'] = query_seq
                feature['processed_ref_seq'] = ref_seq
                processed_features.append(feature)
                
                # 保存比对信息到文件（原始和处理后）
                self._save_feature_alignment(feature)
            except Exception as e:
                print(f"Warning: Failed to process CDS {feature.get('attributes', {}).get('ID')}: {str(e)}")
        
        return processed_features

    def _save_feature_alignment(self, feature):
        """
        保存单个特征的比对结果到文件（包括原始和处理后序列）
        """
        feature_id = feature['attributes'].get('ID', f"feature_{len(os.listdir(self.tmp_dir))}")
        
        # 保存原始比对
        original_align_file = os.path.join(self.tmp_dir, f"{feature_id}_original_alignment.fasta")
        with open(original_align_file, 'w') as f:
            f.write(f">reference\n{feature['original_ref_seq']}\n")
            f.write(f">query\n{feature['original_query_seq']}\n")
        
        # 保存处理后比对
        processed_align_file = os.path.join(self.tmp_dir, f"{feature_id}_processed_alignment.fasta")
        with open(processed_align_file, 'w') as f:
            f.write(f">reference\n{feature['processed_ref_seq']}\n")
            f.write(f">query\n{feature['processed_query_seq']}\n")
        
        # 保存调整信息
        info_file = os.path.join(self.tmp_dir, f"{feature_id}_adjustment_info.txt")
        with open(info_file, 'w') as f:
            original_len = len(feature['original_query_seq'])
            processed_len = len(feature['processed_query_seq'])
            f.write(f"Feature ID: {feature_id}\n")
            f.write(f"Original length: {original_len}\n")
            f.write(f"Processed length: {processed_len}\n")
            f.write(f"Gaps added: {processed_len - original_len}\n")
            f.write(f"Adjustment reason: Ensure length is multiple of 3 for translation\n")

    def translate_cds_features(self, features):
        """
        翻译处理后的CDS序列（负链已处理）
        """
        translations = {}
        
        for feature in features:
            feature_id = feature['attributes'].get('ID', f"cds_{len(translations)}")
            seq = feature['processed_query_seq']
            
            try:
                # 检查起始密码子（所有链方向，因为负链已处理）
                start_codon = str(seq[:3]).upper().replace('-', 'N')
                valid_starts = ['ATG', 'GTG', 'TTG', 'ATA', 'ATT', 'ATC']
                if start_codon not in valid_starts:
                    translations[feature_id] = {
                        'protein_seq': None,
                        'status': f"No valid start codon: {start_codon}",
                        'dna_seq': seq
                    }
                    continue
                
                # 序列方向已正确处理（负链已反向互补）
                to_translate = Seq(str(seq).replace('N', '-').replace('-', 'N'))
                
                # 自定义翻译函数，处理不完整密码子
                def custom_translate(seq):
                    protein = []
                    for i in range(0, len(seq), 3):
                        codon = seq[i:i+3]
                        if len(codon) < 3:
                            # 不完整密码子用?代替
                            protein.append('?')
                            continue
                        try:
                            aa = codon.translate(table=self.genetic_code, cds=False, stop_symbol='*')
                            protein.append(str(aa))
                        except:
                            # 无法翻译的密码子用?代替
                            protein.append('?')
                    return ''.join(protein)
                
                # 翻译（使用指定的遗传密码表）
                protein_seq = custom_translate(to_translate)
                
                # 检查终止密码子
                stop_codon_pos = protein_seq.find('*') if '*' in protein_seq else -1
                if stop_codon_pos >= 0:
                    # 截取到第一个终止密码子
                    protein_seq = protein_seq[:stop_codon_pos+1]
                    status = "Success (truncated at stop codon)"
                    if stop_codon_pos < len(protein_seq) - 1:
                        status += " - WARNING: premature stop codon detected"
                else:
                    status = 'Success (no stop codon found)'
                
                if '?' in protein_seq:
                    status += " (with ambiguous translations)"
                
                translations[feature_id] = {
                    'protein_seq': protein_seq,
                    'status': status,
                    'dna_seq': seq
                }
            except Exception as e:
                translations[feature_id] = {
                    'protein_seq': None,
                    'status': f"Translation failed: {str(e)}",
                    'dna_seq': seq
                }
        
        return translations

    def _write_gff3(self, features):
        """生成GFF3注释文件，正确输出分段CDS"""
        filename = f"{self.prefix}_annotation.gff3" if self.prefix else "annotation.gff3"
        gff_output = os.path.join(self.output_dir, filename)
        
        with open(gff_output, 'w') as f:
            f.write("##gff-version 3\n")
            f.write(f"##sequence-region {self.query_seq.id} 1 {len(self.query_seq.seq)}\n")
            f.write(f"##date {datetime.datetime.now().strftime('%Y-%m-%d')}\n")
            f.write(f"##source MitoAnnotator\n")
            
            # 写入基因组特征
            f.write(f"{self.query_seq.id}\tMitoAnnotator\tregion\t1\t{len(self.query_seq.seq)}\t.\t+\t.\tID={self.query_seq.id};Name=mitochondrion\n")
            
            # 按特征类型排序写入
            feature_order = ['gene', 'mRNA', 'CDS', 'tRNA', 'rRNA', 'control_region']
            for feat_type in feature_order:
                if feat_type in features:
                    for feature in features[feat_type]:
                        # 处理分段CDS
                        if feat_type == 'CDS' and 'segments' in feature:
                            for segment in feature['segments']:
                                attrs = ';'.join(f"{k}={v}" for k, v in segment['attributes'].items())
                                f.write(f"{self.query_seq.id}\tMitoAnnotator\t{feat_type}\t{segment['start']}\t{segment['end']}\t.\t{segment['strand']}\t{segment['phase']}\t{attrs}\n")
                        else:
                            attrs = ';'.join(f"{k}={v}" for k, v in feature['attributes'].items())
                            f.write(f"{self.query_seq.id}\tMitoAnnotator\t{feat_type}\t{feature['start']}\t{feature['end']}\t.\t{feature['strand']}\t{feature.get('phase','.')}\t{attrs}\n")

    def _write_gb(self, features, translations):
        """生成GenBank格式注释文件，处理分段CDS"""
        filename = f"{self.prefix}_annotation.gb" if self.prefix else "annotation.gb"
        gb_output = os.path.join(self.output_dir, filename)
        
        # 创建SeqRecord对象
        record = SeqRecord(
            Seq(str(self.query_seq.seq)),
            id=self.query_seq.id,
            name=self.query_seq.id,
            description=f"Annotated mitochondrial genome",
            annotations={
                "molecule_type": "DNA",
                "topology": "circular",
                "data_file_division": "INV",
                "date": datetime.datetime.now().strftime("%d-%b-%Y").upper(),
                "accessions": [self.query_seq.id],
                "sequence_version": 1,
                "source": "mitochondrion"
            }
        )
        
        # 添加区域特征
        region_feature = SeqFeature.SeqFeature(
            location=SeqFeature.FeatureLocation(0, len(self.query_seq.seq), strand=1),
            type="source",
            qualifiers={
                "organism": self.query_seq.id,
                "organelle": "mitochondrion",
                "mol_type": "genomic DNA"
            }
        )
        record.features.append(region_feature)
        
        # 添加其他特征
        for feat_type, feature_list in features.items():
            for feature in feature_list:
                # 跳过region特征
                if feat_type == "region":
                    continue
                    
                # 处理分段CDS
                if feat_type == "CDS" and 'segments' in feature:
                    # 创建组合位置
                    locations = []
                    for segment in feature['segments']:
                        start = segment['start'] - 1  # 转换为0-based
                        end = segment['end']  # 0-based的end是exclusive，但GenBank使用inclusive end
                        strand = 1 if feature['strand'] == '+' else -1
                        locations.append(SeqFeature.FeatureLocation(start, end, strand=strand))
                    
                    # 创建组合特征位置
                    if len(locations) > 1:
                        # 对于负链，需要反转位置顺序
                        if feature['strand'] == '-':
                            locations.reverse()
                        compound_location = SeqFeature.CompoundLocation(locations)
                    else:
                        compound_location = locations[0]
                    
                    # 创建特征
                    seq_feature = SeqFeature.SeqFeature(
                        location=compound_location,
                        type=feat_type,
                        qualifiers=feature['attributes'].copy()
                    )
                    
                    # 添加翻译产物
                    feature_id = feature['attributes'].get('ID', '')
                    if feature_id in translations and translations[feature_id]['protein_seq']:
                        seq_feature.qualifiers['translation'] = translations[feature_id]['protein_seq']
                        seq_feature.qualifiers['transl_except'] = "contains ambiguous translations" if '?' in translations[feature_id]['protein_seq'] else ""
                        if "premature stop codon" in translations[feature_id]['status']:
                            seq_feature.qualifiers['note'] = "contains premature stop codon"
                    
                    record.features.append(seq_feature)
                else:
                    # 处理非分段特征
                    location = SeqFeature.FeatureLocation(
                        feature['start'] - 1,
                        feature['end'],
                        strand=(1 if feature['strand'] == '+' else -1)
                    )
                    
                    seq_feature = SeqFeature.SeqFeature(
                        location=location,
                        type=feat_type,
                        qualifiers=feature['attributes'].copy()
                    )
                    
                    # 对于CDS特征，添加翻译产物
                    if feat_type == "CDS":
                        feature_id = feature['attributes'].get('ID', '')
                        if feature_id in translations and translations[feature_id]['protein_seq']:
                            seq_feature.qualifiers['translation'] = translations[feature_id]['protein_seq']
                            seq_feature.qualifiers['transl_except'] = "contains ambiguous translations" if '?' in translations[feature_id]['protein_seq'] else ""
                            if "premature stop codon" in translations[feature_id]['status']:
                                seq_feature.qualifiers['note'] = "contains premature stop codon"
                    
                    record.features.append(seq_feature)
        
        # 写入GenBank文件
        with open(gb_output, 'w') as f:
            SeqIO.write(record, f, "genbank")

    def _write_protein_sequences(self, translations):
        """
        输出蛋白质序列
        """
        filename = f"{self.prefix}_proteins.fasta" if self.prefix else "proteins.fasta"
        output_file = os.path.join(self.output_dir, filename)
        
        with open(output_file, 'w') as f:
            for feature_id, data in translations.items():
                if data['protein_seq']:
                    f.write(f">{feature_id}\n")
                    for i in range(0, len(data['protein_seq']), 80):
                        f.write(f"{str(data['protein_seq'][i:i+80])}\n")

    def _write_dna_sequences(self, translations):
        """
        输出CDS的DNA序列
        """
        filename = f"{self.prefix}_cds_sequences.fasta" if self.prefix else "cds_sequences.fasta"
        output_file = os.path.join(self.output_dir, filename)
        
        with open(output_file, 'w') as f:
            for feature_id, data in translations.items():
                seq = str(data['dna_seq']).replace('-', 'N')
                f.write(f">{feature_id}\n")
                for i in range(0, len(seq), 80):
                    f.write(f"{seq[i:i+80]}\n")

    def _generate_report(self, features, translations):
        """
        生成统计报告
        """
        filename = f"{self.prefix}_report.txt" if self.prefix else "report.txt"
        report_output = os.path.join(self.output_dir, filename)
        
        # 计算统计数据
        feature_counts = {k: len(v) for k, v in features.items()}
        translation_stats = {
            'success': sum(1 for t in translations.values() if 'Success' in t['status']),
            'no_start': sum(1 for t in translations.values() if 'No valid start' in t['status']),
            'premature_stop': sum(1 for t in translations.values() if 'premature stop codon' in t['status']),
            'failed': sum(1 for t in translations.values() if t['status'] not in ['Success', 'No valid start'])
        }
        
        with open(report_output, 'w') as f:
            f.write("Mitochondrial Genome Annotation Report\n")
            f.write("=====================================\n\n")
            f.write(f"Sequence ID: {self.query_seq.id}\n")
            f.write(f"Sequence length: {len(self.query_seq.seq)} bp\n")
            f.write(f"Genetic code: {self.genetic_code}\n\n")
            
            f.write("Feature counts:\n")
            for feat_type, count in feature_counts.items():
                f.write(f"  {feat_type}: {count}\n")
            
            f.write("\nTranslation results:\n")
            f.write(f"  Successfully translated: {translation_stats['success']}\n")
            f.write(f"  No valid start codon: {translation_stats['no_start']}\n")
            f.write(f"  Premature stop codons: {translation_stats['premature_stop']}\n")
            f.write(f"  Failed translations: {translation_stats['failed']}\n\n")
            
            f.write("Detailed translation status:\n")
            for feature_id, data in translations.items():
                f.write(f"  {feature_id}: {data['status']}\n")

    def run(self):
        """
        执行整个注释流程
        """
        print(f"Processing {self.query_file}...")
        
        # 检查输入文件
        for f in [self.query_file, self.reference_file, self.reference_annotation]:
            if not os.path.exists(f):
                raise FileNotFoundError(f"Input file not found: {f}")
        
        try:
            # 1. 全序列比对
            ref_aligned, query_aligned = self.run_whole_alignment()
            
            # 2. 提取特征区域
            features = self.extract_aligned_features(ref_aligned, query_aligned)
            
            # 3. 处理CDS特征
            processed_cds = self.process_cds_features(features)
            
            # 4. 翻译CDS
            translations = self.translate_cds_features(processed_cds)
            
            # 5. 生成输出文件
            self._write_gff3(features)
            self._write_gb(features, translations)  # 新增的GenBank输出
            self._write_protein_sequences(translations)
            self._write_dna_sequences(translations)
            self._generate_report(features, translations)
            
            print(f"Successfully annotated {self.query_seq.id}")
            return True
        except Exception as e:
            print(f"Annotation failed: {str(e)}")
            return False
        finally:
            # 保留中间文件（注释掉以下行以删除临时文件）
            # if os.path.exists(self.tmp_dir):
            #     shutil.rmtree(self.tmp_dir)
            pass

def batch_process(query_dir, reference_file, reference_annotation, output_dir, genetic_code=2, prefix=None):
    """
    批量处理目录中的FASTA文件
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取所有FASTA文件
    fasta_files = glob.glob(os.path.join(query_dir, "*.fa")) + \
                 glob.glob(os.path.join(query_dir, "*.fasta")) + \
                 glob.glob(os.path.join(query_dir, "*.fna"))
    
    print(f"Found {len(fasta_files)} FASTA files in {query_dir}")
    
    success_count = 0
    failed_count = 0
    
    for fasta_file in fasta_files:
        file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
        file_output_dir = os.path.join(output_dir, file_prefix)
        os.makedirs(file_output_dir, exist_ok=True)
        
        try:
            annotator = MitochondrialAnnotation(
                query_file=fasta_file,
                reference_file=reference_file,
                reference_annotation=reference_annotation,
                output_dir=file_output_dir,
                genetic_code=genetic_code,
                prefix=file_prefix if prefix is None else f"{prefix}_{file_prefix}"
            )
            if annotator.run():
                success_count += 1
            else:
                failed_count += 1
        except Exception as e:
            print(f"Error processing {fasta_file}: {str(e)}")
            failed_count += 1
    
    print(f"\nBatch processing completed. Success: {success_count}, Failed: {failed_count}")

def main():
    parser = argparse.ArgumentParser(description="Annotate mitochondrial genomes")
    parser.add_argument("--query", required=True, help="Query FASTA file or directory")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--annotation", required=True, help="Reference GFF annotation file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--genetic_code", type=int, default=2, help="Genetic code (default: 2 - Vertebrate Mitochondrial)")
    parser.add_argument("--prefix", help="Prefix for output files")
    parser.add_argument("--batch", action="store_true", help="Batch process files in query directory")
    
    args = parser.parse_args()
    
    if args.batch:
        batch_process(
            query_dir=args.query,
            reference_file=args.reference,
            reference_annotation=args.annotation,
            output_dir=args.output,
            genetic_code=args.genetic_code,
            prefix=args.prefix
        )
    else:
        annotator = MitochondrialAnnotation(
            query_file=args.query,
            reference_file=args.reference,
            reference_annotation=args.annotation,
            output_dir=args.output,
            genetic_code=args.genetic_code,
            prefix=args.prefix
        )
        annotator.run()

if __name__ == "__main__":
    main()